# Hull hydrostatics

The pipeline is `read_bri -> hull -> hydrostatics / solve_float`. It uses base
MATLAB (including `polyshape` and constrained Delaunay triangulation), and is
tested with MATLAB R2025a. Geometry remains independent of loading and pose.
The existing Moth setup and resistance model are not connected automatically.

## Run the independent box example and tests

From the repository root in MATLAB:

```matlab
addpath(pwd)
run('examples/upright_box.m')
results = runtests('tests');
assertSuccess(results)
```

The box is 4 m long, 2 m wide and 1.5 m high. At 0.5 m draft its volume is
4 m^3, buoyancy centre in body coordinates is `[2 0 1.25]` m, waterplane area
is 8 m^2, and wetted area including its bottom, sides and ends is 14 m^2.
Waterplane second moments about its centroid are 8/3 m^4 (roll) and
32/3 m^4 (pitch). At 1000 kg/m^3, its displaced mass is 4000 kg.

## Coordinates and interfaces

Body and world coordinates use x forward, y starboard and z downward. The
world water surface is z=0. Body origins are chosen explicitly at import;
`box_hull` places the origin at the stern centreline on deck.

```matlab
hull = box_hull(4, 2, 1.5);
env = struct('water', struct('rho', 1000), 'g', 9.81);
pose = struct('elevation', 1, 'heel', 0, 'trim', 0);
hs = hydrostatics(hull, pose, env);
```

`pose.elevation` is upward-positive, in metres. It is not draft: for this
upright box, elevation 1 m means draft 0.5 m. Angles are in degrees.
A body column vector transforms as `world = Ry(trim)*Rx(heel)*body - [0;0;elevation]`.
Positive heel lowers the starboard side; positive trim raises the bow.
Rotations pivot about the body origin. Heel is applied first, then trim about
the world y axis; the instantaneous heel axis after trim is `Ry*[1;0;0]`.

A manually supplied hull requires `vertices` (N-by-3 metres) and `faces`
(M-by-3 one-based triangle indices). Faces must form one connected, closed,
consistently wound shell. Global reversal is accepted and corrected locally;
mixed winding, unused vertices, degenerate triangles and open edges are errors.
The input struct is not mutated. Imported hulls also contain `name`, `units`,
`axes`, `sections` and `source` metadata.

`hydrostatics` returns:

| Field | Meaning / units |
| --- | --- |
| `volume`, `displacedMass`, `buoyancyForce` | m^3, kg, upward force magnitude in N |
| `centreOfBuoyancyBody`, `centreOfBuoyancyWorld` | Volume centroid, m |
| `wettedArea` | Submerged original hull surface, m^2; includes sealed end caps |
| `waterplaneArea`, `centreOfFlotationWorld` | Waterplane area (m^2) and centroid (m) |
| `waterplaneMoments.roll`, `.pitch`, `.product` | Centroidal integrals of y^2, x^2, xy, m^4 |
| `draft` | Maximum positive world depth of the hull, m |
| `waterlineExtent`, `LWL` | Min/max world x at the waterplane and their difference, m |
| `condition` | `dry`, `partially_submerged`, or `fully_submerged` |
| `mesh.vertices`, `.faces`, `.isWaterplane` | Clipped surface and artificial cap diagnostics |

`LWL` is a world-x projected extent, not length along a tilted body axis.
Dry hulls return zero displaced volume and NaN buoyancy centroids. No-waterplane
conditions return zero area/moments and NaN horizontal flotation coordinates.
Exact bottom tangency is classified dry; exact top tangency is fully submerged.
At top tangency the full closed shell area is counted, and waterplane quantities
are zero. Use a slightly immersed/exposed pose for one-sided waterplane limits.

Clipped diagnostic vertices can be duplicated at cap boundaries; the mesh is
geometrically closed but is not intended as a welded export mesh. Waterplane
caps carry no wetted area. Volume/centroid use signed tetrahedral integration;
waterplane moments use triangle integration. No perimeter-times-length estimate
is used for wetted area.

## Import the section-block BRI dialect

```matlab
options = struct( ...
    'lengthScale', 1, ...           % source lengths -> metres
    'xHeaderColumn', 2, ...         % 2 or 3 in [count value value]
    'axisSigns', [1 1 -1], ...      % source xyz sign conversion
    'originOffset', [0 0 0], ...    % metres, added AFTER scale/sign conversion
    'halfHull', true);
% Use only after deciding these settings for the source file.
% The supplied kayak also requires an explicit sealed approximation:
options.halfHullClosure = 'seal';
hull = read_bri('path/to/verified_geometry.bri', options);
```

The five mapping options are mandatory: the example values above are provisional,
not verified kayak conventions. Import assumes coordinate pairs are y,z; arbitrary
axis permutations are unsupported. Lengths map by
`body = source .* (lengthScale*axisSigns) + originOffset`.

Supported file layout: name line, one numeric file-header value, repeated
`pointCount headerValue2 headerValue3` blocks containing y/z pairs and a numeric
scalar trailer, followed by `0 0 0`. Blank lines and arbitrary whitespace are
accepted; non-finite numbers, inconsistent counts, missing terminators and trailing
data are rejected. Header values and trailers are preserved without interpreting
flags. Descending stations are reversed; unordered or repeated stations are errors.

Half-sections must stay on one side of y=0. Default `halfHullClosure='strict'`
requires both endpoints on y=0. Explicit `halfHullClosure='seal'` permits open
endpoints, connecting each to its mirrored partner by a straight transverse edge,
and removes zero-area centreline tails. This creates a sealed approximation,
including an artificial deck over any cockpit opening; it is not downflooding
geometry. Original pairs remain in `sections(k).sourcePoints`, with changes
flagged by `closureAdjusted`. The supplied kayak has open endpoints and also
returns internally to the centreline at station 1.4 before extending outward
again. Under these provisional settings, even sealed import correctly rejects
that station as an intersecting/touching contour (edges 75 and 108). Run
`run('examples/kayak_smoke.m')` to reproduce the diagnostic. No physical kayak
hydrostatics are validated; corrected sections or confirmed export semantics
are needed. Half-sections are mirrored before coordinate conversion.
Consecutive duplicate points are removed;
self-intersecting or zero-area contours are rejected. Fully closed nondegenerate
sections are required: point/line tips must first be represented by suitable
nondegenerate sections. There is no automatic bow/stern extrapolation.

Each section is oriented consistently and its seam is chosen at minimum y,
breaking ties by minimum z. A shared normalized-perimeter grid retains the union
of section vertex fractions (merged within 1e-12). Adjacent rings are triangulated;
ends are sealed with constrained planar triangulation. This is a piecewise planar
loft, not a fitted smooth surface. Inspect correspondence on shapes with abrupt
changes of section form. Input mesh validation checks topology and winding, not
arbitrary three-dimensional triangle self-intersections.

## Floating equilibrium

```matlab
loading = struct('mass', 4000, 'cg', [2 0 1]);
upright = solve_float(hull, loading, struct('mode','draft'), env);
heeled = solve_float(hull, loading, ...
    struct('mode','draft','heel',10,'trim',0), env);
free = solve_float(hull, loading, ...
    struct('mode','equilibrium','heel',3,'trim',2), env);
assert(free.converged, free.message)
```

Use a box `hull` with the illustrated box loading. Do not combine the kayak
geometry with Moth or box loading by implication.

`mode='draft'` (default) solves vertical position at fixed angles using bracketed
`fzero`. Bounds come from rotated hull extents. Nonpositive mass and mass at or
above maximum sealed displacement are rejected. Successful solutions require
relative mass residual below 1e-8. `loading.cg` is required in both modes to report
moment residuals; fixed-attitude draft convergence does not imply moment balance.

`mode='equilibrium'` also aligns horizontal buoyancy and gravity centres using a
finite-difference damped Newton solve. Defaults: initial angles zero, maximum
30 iterations, bounds `[-30 30; -30 30]` (heel row, trim row). Bounds may be changed
but must stay strictly inside +/-89 degrees. Steps outside bounds are rejected.
The finite-difference increment is 0.001 degrees; horizontal convergence tolerance
is 1e-8 times hull scale. Root position tolerance is 1e-11 times hull scale.

Results include `pose`, `hydrostatics`, `converged`, `message`, `iterations`,
`massResidual` (kg), `forceResidualWorld` (N), `horizontalResidual` (m), and
`momentResidualWorld` (N m, buoyancy moment about CG). Non-convergence returns
`converged=false` and the last evaluated state; callers must inspect that flag.
Malformed geometry and failed inner draft tolerances raise errors.

At converged free equilibrium, `restoringStiffness` is the symmetrized negative
Jacobian of generalized roll/pitch moments in N m/rad, with draft re-solved at
every perturbation. Its eigenvalues classify `stable`, `unstable`, or
`near_neutral`, using threshold `1e-6*mass*g*hullScale`. Small diagnostic
perturbations may extend beyond solver angle bounds. Other states report
`not_evaluated`. This is a local classification, not proof of global uniqueness.

## Numerical and physical limits

Clipping snaps depths within 1e-10 times hull scale to the water surface and
welds nearby clipped vertices using the same tolerance. Exact waterplane contacts
creating branching contours are rejected with a request to perturb the pose.
Very thin features or extremely small displacements may require a better-scaled
model. Refinement tests use a circular cylinder and confirm decreasing errors
against its analytic half-submerged volume and wetted area.

The model is a sealed, rigid hull in calm water. Flooding, downflooding, tank
free surfaces, foil lift, sails and external loads are not included. Mesh
refinement and source-coordinate verification remain necessary for real hulls.
