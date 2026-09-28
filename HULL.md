# Hull resistance

The standalone `delft_resistance` module combines the **Keuning–Katgert (2008)**
upright, untrimmed bare-hull regression (equation 1.7, Table 2) with ITTC-1957
friction using waterline length. It accepts fixed, static immersed geometry.
It does not calculate hydrostatics, heel, leeway, appendage drag, foil support,
sail-induced trim, roughness or added resistance in waves. “Untrimmed” here
means no applied sail-induced trimming moment; the regression incorporates the
behaviour of the tested hulls at speed, rather than solving dynamic trim itself.

## Calling interface

```matlab
[R, details] = delft_resistance(V, hull, env);
```

`V` is speed **through water**, in m/s: a nonempty, real, nonnegative numeric
scalar, row vector or column vector. Numeric results are double arrays with
the same shape. `R` is the resistance scalar for each speed in newtons, not a
body-axis force vector. For a forward-moving hull, the corresponding drag
force acts in negative x. The regression is not clamped: inspect the negative
result flags before treating its output as a physically meaningful magnitude.

Required hull fields are finite real numeric scalars:

| Field | Meaning | Unit |
| --- | --- | --- |
| `LWL` | Static waterline length | m |
| `BWL` | Maximum waterline beam | m |
| `draft` | Canoe-body draft, excluding appendages | m |
| `volume` | Canoe-body displacement volume, not mass | m³ |
| `Swet` | Bare-hull wetted surface area | m² |
| `Awp` | Static waterplane area | m² |
| `Cp` | Prismatic coefficient | dimensionless |
| `Cm` | Midship section coefficient | dimensionless |
| `xLCB` | Longitudinal centre of buoyancy | m |
| `xLCF` | Longitudinal centre of flotation | m |
| `xFPP` | Forward endpoint of the static waterline | m |

The environmental structure matches `get_env_params`: `env.water.rho`
(kg/m³), `env.water.nu` (m²/s) and `env.g` (m/s²), all finite and positive.
The caller supplies mutually consistent hull properties for one loading
condition. The module checks basic physical and regression bounds; it cannot
verify the supplied hydrostatics without hull geometry.

### Coordinate convention

Use a right-handed frame with origin at the stern on deck, x forward,
y starboard and z downward. All three longitudinal coordinates must use the
same origin. `draft` is an immersed depth, not a z coordinate from the deck.
The aft waterline endpoint is `xFPP - LWL`; it need not coincide with the origin.
Only longitudinal differences enter this module:

```matlab
LCB_fpp = hull.xFPP - hull.xLCB;
LCF_fpp = hull.xFPP - hull.xLCF;
```

These are positive distances aft from the forward waterline endpoint, as
required by the source equation. Translating all three coordinates by the
same amount must not change resistance.

## Calculation

Using `L = hull.LWL`, `S = hull.Swet`, and `volume = hull.volume`:

```text
Fn = V / sqrt(g * L)
Re = V * L / nu
Cf = 0.075 / (log10(Re) - 2)^2
Rf = 0.5 * rho * V^2 * S * Cf

q = a1 * LCB_fpp/L + a2 * Cp + a3 * volume^(2/3)/Awp
  + a4 * BWL/L + a5 * LCB_fpp/LCF_fpp + a6 * BWL/draft + a7 * Cm

Rr = rho * g * volume * (a0 + q * volume^(1/3)/L)
R  = Rf + Rr
```

The constant `a0` is outside the final length-ratio multiplier. Table 2's
coefficients need no additional factor of 1000. The private coefficient table
contains 13 rows at `Fn = 0.15:0.05:0.75`. Each coefficient is interpolated
linearly between rows. There is no extrapolation or smoothing spline.

Friction retains this repository's full-`LWL` Reynolds length, with no form
factor. This choice is explicit: the combined total should not be described as
an exact reproduction of every historical Delft friction-extrapolation
convention. A fully turbulent wetted surface is assumed. The `Re > 100` guard
only excludes the singularity and lower mathematical branch of the formula;
it is **not** a turbulence or experimental-validity criterion.

## Applicability and diagnostics

Nonzero speed must satisfy `0.15 <= Fn <= 0.75`, including endpoints. The
interval `0 < Fn < 0.15` is unsupported. A few double-precision ulps at either
endpoint are tolerated only to avoid rejecting arithmetic roundoff.

At `V = 0`, `R`, `Rf`, `Rr`, `Fn` and `Re` are zero; `Cf` is `NaN` because it
is not evaluated. Hull and environment validation still apply.

The following geometry screen uses the **printed marginal extrema in the
2008 Table 1**, not the older 1998 limits:

| Ratio or coefficient | Minimum | Maximum |
| --- | ---: | ---: |
| `LCB_fpp/LWL` | 0.500 | 0.582 |
| `Cp` | 0.519 | 0.599 |
| `volume^(2/3)/Awp` | 0.079 | 0.265 |
| `BWL/LWL` | 0.170 | 0.366 |
| `LCB_fpp/LCF_fpp` | 0.920 | 1.002 |
| `volume^(1/3)/LWL` | 0.12 | 0.23 |
| `Cm` | 0.646 | 0.790 |
| `BWL/draft` | 2.46 | 19.38 |

These rounded extrema are an implementation screen, **not a validated joint
hull-design envelope**. Passing them does not establish that a combination of
parameters was tested. More precise values near the rounded limits can be
rejected. The paper also uses fewer models above `Fn = 0.60`; this module does
not claim speed-specific geometric validation.

`details` contains:

- `Rf`, `Rr`, `Fn`, `Re`, `Cf`: arrays shaped like `V`.
- `model`: identifies the regression and friction convention.
- `validity.withinSpeedRange`: true for evaluated nonzero speeds; false at zero.
- `validity.zeroSpeed`: marks the explicit zero-speed case.
- `validity.geometry`: screen basis, ratio names, input values and limits.
- `validity.negativeResiduary`, `validity.negativeTotal`: flag unmodified negative results.
- `validity.reducedHighSpeedDataset`: marks speeds above `Fn = 0.60`.
- `validity.note`: states the limits of the applicability screen.

A negative fitted residual is preserved for diagnosis, not interpreted as
thrust. Missing fields, malformed inputs, impossible geometry, unsupported
ratios, unsupported speed or nonfinite calculations raise errors. An invalid
speed rejects the entire call, with its linear index in the message.
Error identifiers are `delft:MissingField`, `delft:InvalidInput`,
`delft:InvalidSpeed`, `delft:InvalidGeometry`, `delft:GeometryRange`,
`delft:SpeedRange`, `delft:ReynoldsRange` and `delft:NumericalRange`.

## Reference example and verification

From the repository root:

```matlab
addpath('examples');
[hull, env] = delft_reference_hull();
V = [0, 0.15, 0.425, 0.60, 0.75] * sqrt(env.g * hull.LWL);
[R, details] = delft_resistance(V, hull, env);
example_delft_resistance(); % plot total, friction and residuary components

results = runtests('tests');
assertSuccess(results);
```

The example reconstructs **DSYHS hull 25 at LWL = 10 m**. The dimensions,
volume and wetted area in de Baar et al., Table 2, are scaled from 2 m by a
linear factor of 5 (areas by 25, volume by 125). `Cp`, `Cm`, LCB/LWL,
LCB/LCF and volume^(2/3)/Awp come from Keuning–Katgert Table 1. Their rounding
limits reconstruction accuracy. No Moth parameters are substituted or changed.
The example's `xFPP = 12 m` is an illustrative coordinate offset, not a
published deck dimension. Water density 1025 kg/m³, viscosity 1.19e-6 m²/s and
gravity 9.81 m/s² are explicit example assumptions.

Independent 40-digit scalar arithmetic supplies regression references at all
13 coefficient nodes and an interior point. At `Fn = 0.425`, expected
`Rr = 1135.456548478192 N` and `Rf = 371.095530156165 N`. Numerical test
tolerances are 1e-8 N for residuary resistance and 1e-9 N for friction.
These tight tolerances check implementation arithmetic, not physical accuracy.

For a separate graphical comparison, the measured square in **2008 Figure 6**
at `Fn = 0.60` is approximately **3500 N residuary resistance**. The module
returns **3490.879 N** for the reconstruction. The test uses a **±150 N**
comparison allowance: the figure has 500 N ticks, its marker is read visually,
the hull ratios are rounded, and water properties are assumed. This allowance
is not a quantified experimental uncertainty or a general model error bound.
No claim of total-resistance validation follows from this residuary comparison.

Tests also cover vector/scalar equivalence, array shape, summation, coordinate
translation, dimensional scaling, density scaling, zero speed, range boundaries,
negative residual preservation, invalid inputs and the existing environment
structure. The example and tests require base MATLAB only.

## Load a BRI hull into resistance

From the repository root in MATLAB:

```matlab
addpath(pwd)
addpath('examples')
result = kayak_resistance(130); % total kayak + paddler + gear mass, kg
```

The example performs `read_bri -> solve_float -> resistance_hull ->
delft_resistance`. It returns the extracted properties and the rejection
message if the hull fails the Delft geometry screen. It never extrapolates.
To use a different file or confirmed import settings, the underlying calls are:

```matlab
mesh = read_bri(filename, importOptions); % explicit, verified source mapping
env = get_env_params();
loading = struct('mass', totalMass, 'cg', cgBody); % kg and body-frame metres
floating = solve_float(mesh, loading, ...
    struct('mode','draft','heel',0,'trim',0), env);
assert(floating.converged, floating.message)
[resHull, hs] = resistance_hull(mesh, floating.pose, env);
V = [0 .15 .25 .35 .45] * sqrt(env.g * resHull.LWL);
[R, details] = delft_resistance(V, resHull, env);
```

`filename`, `importOptions`, `totalMass` and `cgBody` are user inputs.
Fixed-attitude draft solving balances mass only. The example uses a clearly
labelled placeholder CG, which does not affect its solved draft; it does not
establish moment balance. The adapter rejects nonzero heel/trim and dry or
fully submerged states. A prescribed upright pose can also be passed directly
to `resistance_hull` when the waterline is known.

The adapter measures maximum waterline beam from the waterplane cap and
intersects the immersed triangular mesh at the midpoint of the waterline.
With that immersed midship area `Am`, it uses `Cp=volume/(LWL*Am)` and
`Cm=Am/(BWL*draft)`. This convention is explicit: for unusual shapes where
the reference/maximal section differs from midship, review the coefficient
definition before using a regression. It does not silently substitute a
sampled maximum section. All longitudinal coordinates retain the mesh origin.

### Provisional 130 kg kayak check (2026-09-28)

MATLAB R2025a, repository freshwater density 997.8 kg/m³, zero heel/trim:

| Quantity | Calculated value |
| --- | ---: |
| Draft | 0.128583 m |
| Displacement volume | 0.130287 m³ |
| Waterline length | 2.9600 m |
| Waterline beam | 0.6026 m |
| Wetted surface | 1.7246 m² |
| Waterplane area | 1.4145 m² |
| Cp (midship convention) | 0.7213 |
| Cm | 0.7875 |
| LCB distance aft of FPP / LWL | 0.460953 |

The first rejection is `LCB_fpp/LWL < 0.500`. Independently, `Cp > 0.599`.
No total resistance curve is returned for this kayak. Reversing x alone
would not fix the excessive Cp. Confirm the source mapping and shape before
deciding whether another resistance model is needed.

These are numerical results for the provisional sealed mesh, not validated
kayak predictions: the source mapping remains unconfirmed, all 16 sections
require closure adjustments, and section correspondence affects the loft.
The modelled waterline spans the complete imported station length. Inspect
the shape and measured dimensions/waterline before physical use.

`TestResistanceIntegration` checks analytical box properties, a tapered
homothetic trapezoidal hull with analytical volume/waterplane/coefficients,
successful resistance evaluation, mass-to-draft coupling, coordinate
translation, unsupported poses and the actual 130 kg BRI rejection path.
The synthetic hull is a numerical fixture, not physical validation of Delft.

## Sources

1. Keuning, J.A. & Katgert, M. (2008). *A bare hull resistance prediction method
   derived from the results of the Delft Systematic Yacht Hull Series extended
   to higher speeds*. INNOVSAIL, pp. 13–21. Equation 1.7 and Table 2 (PDF p.6),
   Table 1 (PDF p.5), Figure 6 (PDF p.7).
   [TU Delft record](https://research.tudelft.nl/en/publications/a-bare-hull-resistance-prediction-method-derived-from-the-results/),
   [original paper hosted on Scribd](https://www.scribd.com/document/346520042/Keuning-2008-pdf).
   Equation, coefficient signs and plot were visually inspected during implementation.
2. de Baar, J., Roberts, S., Dwight, R. & Mallol, B. (2015).
   *Uncertainty quantification for a sailing yacht hull, using multi-fidelity
   kriging*. Computers & Fluids 123, 185–201.
   [DOI](https://doi.org/10.1016/j.compfluid.2015.10.004),
   [accepted manuscript, Table 2, p.28](https://research.tudelft.nl/files/28596624/main_dsyhs_revision.pdf).

The coefficient data are stored locally; calculations and tests need no network
access. Full source papers are not redistributed in this repository.
