function tests = test_delft_resistance
% Numerical references independently evaluated with 40-digit Decimal
% arithmetic from Table 2 coefficient rows (not the production helper).
tests = functiontests(localfunctions);
end

function setupOnce(t)
root = fileparts(fileparts(mfilename('fullpath')));
t.applyFixture(matlab.unittest.fixtures.PathFixture(root));
t.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(root, 'examples')));
[t.TestData.h, t.TestData.e] = delft_reference_hull();
end

function testAllCoefficientNodes(t)
h = t.TestData.h; e = t.TestData.e;
fn = (15:5:75)/100;
[~, d] = delft_resistance(fn*sqrt(e.g*h.LWL), h, e);
expected = [0.111120285304 13.330238509751 49.587027672389 ...
    123.896612258572 269.964564569685 731.680583469780 ...
    1539.232513486604 2417.558187404255 3035.546134100488 ...
    3490.878543329057 3801.282138786687 4118.792840314191 4284.773200450765];
verifyEqual(t, d.Rr, expected, 'AbsTol', 1e-8);
end

function testIndependentInteriorAndFriction(t)
h = t.TestData.h; e = t.TestData.e;
[R, d] = delft_resistance(.425*sqrt(e.g*h.LWL), h, e);
verifyEqual(t, d.Rr, 1135.45654847819175, 'AbsTol', 1e-8);
verifyEqual(t, d.Rf, 371.09553015616508, 'AbsTol', 1e-9);
verifyEqual(t, R, 1506.55207863435683, 'AbsTol', 1e-8);
end

function testPublishedFigure6(t)
% Keuning-Katgert 2008, Fig.6 (Sysser25): at Fn=.60 the measured square
% is approximately 3500 N. +/-150 N is a graphical comparison allowance
% (500 N tick spacing, rounded geometry, assumed water), NOT a precision
% tolerance or a quantified experimental uncertainty. See HULL.md.
h = t.TestData.h; e = t.TestData.e;
[~, d] = delft_resistance(.60*sqrt(e.g*h.LWL), h, e);
verifyEqual(t, d.Rr, 3500, 'AbsTol', 150);
end

function testScalarAndVectorShapes(t)
h = t.TestData.h; e = t.TestData.e;
v = [0 .15 .425 .75]*sqrt(e.g*h.LWL);
[r, d] = delft_resistance(v, h, e);
[c, dc] = delft_resistance(v', h, e);
verifySize(t, r, [1 4]); verifySize(t, c, [4 1]);
verifyEqual(t, c, r');
for f = {'Rf','Rr','Fn','Re','Cf'}
    verifyEqual(t, dc.(f{1}), d.(f{1})');
end
for k = 1:numel(v)
    verifyEqual(t, delft_resistance(v(k),h,e), r(k), 'AbsTol', 1e-9);
end
verifyEqual(t, r, d.Rf+d.Rr);
verifyTrue(t, d.validity.geometry.passed);
verifyEqual(t, d.validity.zeroSpeed, [true false false false]);
verifyEqual(t, d.validity.reducedHighSpeedDataset, [false false false true]);
end

function testZeroSpeed(t)
[r,d] = delft_resistance([0;0], t.TestData.h, t.TestData.e);
verifyEqual(t, r, [0;0]); verifyEqual(t,d.Rf,[0;0]);
verifyEqual(t,d.Rr,[0;0]); verifyEqual(t,d.Re,[0;0]);
verifyTrue(t,all(isnan(d.Cf)));
verifyFalse(t,any(d.validity.withinSpeedRange));
end

function testCoordinateTranslation(t)
h = t.TestData.h; e = t.TestData.e;
v = [.2 .4 .7]*sqrt(e.g*h.LWL);
r = delft_resistance(v,h,e);
for name = {'xLCB','xLCF','xFPP'}
    h.(name{1}) = h.(name{1}) - 25;
end
verifyEqual(t,delft_resistance(v,h,e),r,'AbsTol',1e-8);
end

function testScalingAndDensity(t)
h = t.TestData.h; e = t.TestData.e;
v = .5*sqrt(e.g*h.LWL);
[r,d] = delft_resistance(v,h,e);
e.water.rho = 2*e.water.rho;
verifyEqual(t,delft_resistance(v,h,e),2*r,'AbsTol',1e-8);
e = t.TestData.e;
for name = {'LWL','BWL','draft','xLCB','xLCF','xFPP'}
    h.(name{1}) = 2*h.(name{1});
end
h.volume = 8*h.volume; h.Swet = 4*h.Swet; h.Awp = 4*h.Awp;
[~, scaled] = delft_resistance(sqrt(2)*v,h,e);
verifyEqual(t,scaled.Rr,8*d.Rr,'AbsTol',1e-8);
end

function testStrictSpeedRange(t)
h = t.TestData.h; e = t.TestData.e; scale = sqrt(e.g*h.LWL);
verifyTrue(t,all(isfinite(delft_resistance([.15 .75]*scale,h,e))));
for fn = [.149999 .750001 .01]
    verifyError(t,@()delft_resistance([.4 fn]*scale,h,e),'delft:SpeedRange');
end
try
    delft_resistance([.4 .1]*scale,h,e);
catch err
    verifySubstring(t,err.message,'V(2)');
end
end

function testInvalidSpeeds(t)
h = t.TestData.h; e = t.TestData.e;
for bad = {[], ones(2), -1, NaN, Inf, 1+1i, '2', true}
    verifyError(t,@()delft_resistance(bad{1},h,e),'delft:InvalidSpeed');
end
end

function testMissingAndMalformedHull(t)
h = t.TestData.h; e = t.TestData.e;
for name = fieldnames(h)'
    verifyError(t,@()delft_resistance(0,rmfield(h,name{1}),e),'delft:MissingField');
end
for bad = {NaN, Inf, [1 2], 1i, '3', true}
    h = t.TestData.h; h.Swet = bad{1};
    verifyError(t,@()delft_resistance(0,h,e),'delft:InvalidInput');
end
h = t.TestData.h; h.draft = 0;
verifyError(t,@()delft_resistance(0,h,e),'delft:InvalidInput');
verifyError(t,@()delft_resistance(0,[],e),'delft:InvalidInput');
end

function testInvalidEnvironment(t)
h = t.TestData.h; e = t.TestData.e;
verifyError(t,@()delft_resistance(0, h, rmfield(e,'g')),'delft:MissingField');
verifyError(t,@()delft_resistance(0, h, struct()),'delft:MissingField');
e.water.nu = -1;
verifyError(t,@()delft_resistance(0,h,e),'delft:InvalidInput');
e = t.TestData.e; e.water.rho = NaN;
verifyError(t,@()delft_resistance(0,h,e),'delft:InvalidInput');
e = t.TestData.e; e.water = [];
verifyError(t,@()delft_resistance(0,h,e),'delft:InvalidInput');
end

function testImpossibleGeometry(t)
h = t.TestData.h; e = t.TestData.e;
h.xLCF = h.xFPP;
verifyError(t,@()delft_resistance(0,h,e),'delft:InvalidGeometry');
h = t.TestData.h; h.Awp = 2*h.LWL*h.BWL;
verifyError(t,@()delft_resistance(0,h,e),'delft:InvalidGeometry');
end

function testGeometryScreen(t)
h = t.TestData.h; e = t.TestData.e;
h.Cp = .7;
verifyError(t,@()delft_resistance(0,h,e),'delft:GeometryRange');
h = t.TestData.h; h.Cm = .8;
verifyError(t,@()delft_resistance(0,h,e),'delft:GeometryRange');
h = t.TestData.h; h.Cp = .519;
verifyEqual(t,delft_resistance(0,h,e),0);
h.Cp = .599;
verifyEqual(t,delft_resistance(0,h,e),0);
end

function testNegativeResidualNotClipped(t)
% Valid marginal inputs can still produce negative fitted residuals.
h = t.TestData.h; e = t.TestData.e; h.Cp = .599;
[~,d] = delft_resistance(.15*sqrt(e.g*h.LWL),h,e);
verifyLessThan(t,d.Rr,0);
verifyTrue(t,d.validity.negativeResiduary);
end

function testReynoldsGuard(t)
h = t.TestData.h; e = t.TestData.e;
e.water.nu = 10;
verifyError(t,@()delft_resistance(.4*sqrt(e.g*h.LWL),h,e),'delft:ReynoldsRange');
end

function testExistingEnvironment(t)
[r,d] = delft_resistance(4,t.TestData.h,get_env_params());
verifyTrue(t,isfinite(r)); verifyGreaterThan(t,d.Rf,0);
end
