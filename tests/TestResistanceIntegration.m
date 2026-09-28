classdef TestResistanceIntegration < matlab.unittest.TestCase
    methods(TestClassSetup)
        function projectPath(t)
            t.applyFixture(matlab.unittest.fixtures.PathFixture(fileparts(fileparts(mfilename('fullpath')))));
        end
    end
    methods(Test)
        function analyticBoxAndRejection(t)
            env=get_env_params(); p=struct('heel',0,'trim',0,'elevation',1);
            h=resistance_hull(box_hull(),p,env);
            t.verifyEqual([h.LWL h.BWL h.draft h.volume h.Swet h.Awp h.Cp h.Cm], ...
                [4 2 .5 4 14 8 1 1],'AbsTol',1e-9);
            t.verifyEqual([h.xLCB h.xLCF h.xFPP h.midshipArea],[2 2 4 1],'AbsTol',1e-9);
            t.verifyError(@()delft_resistance(1,h,env),'delft:GeometryRange');
        end
        function syntheticHullThroughBothModules(t)
            % Analytical tapered trapezoidal hull, constructed only as a
            % numerical integration fixture, not an experimental reference.
            L=10; B=2.5; T=.7; q=[.45 1 .45];
            for k=1:3
                sections(k).x=(k-1)*L/2; %#ok<AGROW>
                sections(k).yz=[-B/2 -1;B/2 -1;B/2 0;.2*B T;-.2*B T;-B/2 0]*q(k); %#ok<AGROW>
            end
            [mesh.vertices,mesh.faces]=hgeom.section_mesh(sections);
            env=get_env_params(); p=struct('heel',0,'trim',0,'elevation',0);
            [h,hs]=resistance_hull(mesh,p,env);
            cp=(1+.45+.45^2)/3;
            t.verifyEqual(h.Cp,cp,'AbsTol',1e-9);
            t.verifyEqual(h.Cm,.7,'AbsTol',1e-9);
            t.verifyEqual(h.volume,cp*L*.7*B*T,'AbsTol',1e-9);
            t.verifyEqual(h.Awp,.725*L*B,'AbsTol',1e-9);
            V=[0 .15 .4 .75]*sqrt(env.g*L);
            [R,d]=delft_resistance(V,h,env);
            t.verifyTrue(all(isfinite(R))); t.verifyEqual(R,d.Rf+d.Rr);
            expected=.5*env.water.rho*h.Swet*V(3)^2*.075/(log10(V(3)*L/env.water.nu)-2)^2;
            t.verifyEqual(d.Rf(3),expected,'RelTol',1e-12);
            load=struct('mass',hs.displacedMass,'cg',[5 0 0]);
            floating=solve_float(mesh,load,struct('mode','draft'),env);
            t.verifyTrue(floating.converged);
            loaded=resistance_hull(mesh,floating.pose,env);
            t.verifyEqual(delft_resistance(V,loaded,env),R,'AbsTol',1e-6);
            mesh.vertices=mesh.vertices+[17 0 3]; p.elevation=3;
            shifted=resistance_hull(mesh,p,env);
            t.verifyEqual(delft_resistance(V,shifted,env),R,'AbsTol',1e-8);
        end
        function unsupportedStates(t)
            env=get_env_params(); mesh=box_hull();
            p=struct('heel',1,'trim',0,'elevation',1);
            t.verifyError(@()resistance_hull(mesh,p,env),'hull:ResistancePose');
            p.heel=0; p.trim=1;
            t.verifyError(@()resistance_hull(mesh,p,env),'hull:ResistancePose');
            p.trim=0; p.elevation=2;
            t.verifyError(@()resistance_hull(mesh,p,env),'hull:ResistanceWaterplane');
            p.elevation=-1;
            t.verifyError(@()resistance_hull(mesh,p,env),'hull:ResistanceWaterplane');
        end
        function exploratoryKayak(t)
            root=fileparts(fileparts(mfilename('fullpath')));
            t.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(root,'examples')));
            t.applyFixture(matlab.unittest.fixtures.SuppressedWarningsFixture('delft:GeometryExtrapolation'));
            previous=get(groot,'DefaultFigureVisible');
            cleanup=onCleanup(@()set(groot,'DefaultFigureVisible',previous)); %#ok<NASGU>
            set(groot,'DefaultFigureVisible','off');
            r=kayak_resistance(130,struct('geometryPolicy','exploratory'));
            figCleanup=onCleanup(@()close(r.figure)); %#ok<NASGU>
            t.verifyTrue(all(isfinite(r.resistance)));
            t.verifyEmpty(r.delftError);
            t.verifyTrue(r.details.validity.geometryExtrapolation);
            g=r.details.validity.geometry;
            t.verifyEqual(g.withinRange,g.values>=g.lower & g.values<=g.upper);
            t.verifyFalse(g.withinRange(1)); t.verifyFalse(g.withinRange(2));
            t.verifyEqual(r.resistance,r.details.Rf+r.details.Rr);
        end
        function sectionAreas(t)
            env=get_env_params(); p=struct('heel',0,'trim',0,'elevation',1);
            hs=hydrostatics(box_hull(),p,env);
            for x=[.01 1 2 3 3.99]
                t.verifyEqual(hgeom.immersed_section_area(hs.mesh,x),1,'AbsTol',1e-9);
            end
        end
        function kayakLoadedPipeline(t)
            root=fileparts(fileparts(mfilename('fullpath')));
            t.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(root,'examples')));
            r=kayak_resistance(130);
            t.verifyTrue(r.floating.converged);
            t.verifyEqual(r.hydrostatics.displacedMass,130,'RelTol',1e-8);
            t.verifyGreaterThan(r.hull.midshipArea,0);
            t.verifyEmpty(r.resistance);
            t.verifyNotEmpty(r.delftError);
            t.verifyError(@()delft_resistance(r.speed,r.hull,get_env_params()),'delft:GeometryRange');
        end
    end
end
