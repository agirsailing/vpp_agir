classdef TestHydrostatics < matlab.unittest.TestCase
    properties
        Hull
        Env
        Pose
    end
    methods(TestClassSetup)
        function projectPath(t)
            t.applyFixture(matlab.unittest.fixtures.PathFixture(fileparts(fileparts(mfilename('fullpath')))));
        end
    end
    methods(TestMethodSetup)
        function setup(t)
            t.Hull=box_hull();
            t.Env=struct('water',struct('rho',1000),'g',9.81);
            t.Pose=struct('elevation',1,'heel',0,'trim',0);
        end
    end
    methods(Test)
        function exactBox(t)
            h=hydrostatics(t.Hull,t.Pose,t.Env);
            t.verifyEqual(h.volume,4,'RelTol',1e-9);
            t.verifyEqual(h.centreOfBuoyancyBody,[2 0 1.25],'AbsTol',1e-9);
            t.verifyEqual(h.waterplaneArea,8,'RelTol',1e-9);
            t.verifyEqual(h.wettedArea,14,'RelTol',1e-9);
            t.verifyEqual(h.waterplaneMoments.roll,8/3,'RelTol',1e-9);
            t.verifyEqual(h.waterplaneMoments.pitch,32/3,'RelTol',1e-9);
            t.verifyEqual(h.waterplaneMoments.product,0,'AbsTol',1e-9);
            t.verifyEqual(h.LWL,4,'AbsTol',1e-9);
            t.verifyEqual(h.buoyancyForce,39240,'RelTol',1e-9);
        end
        function dryAndFull(t)
            p=t.Pose; p.elevation=2; h=hydrostatics(t.Hull,p,t.Env);
            t.verifyEqual(h.condition,'dry'); t.verifyEqual(h.volume,0);
            t.verifyTrue(all(isnan(h.centreOfBuoyancyBody)));
            p.elevation=-1; h=hydrostatics(t.Hull,p,t.Env);
            t.verifyEqual(h.condition,'fully_submerged');
            t.verifyEqual(h.volume,12,'AbsTol',1e-9);
            t.verifyEqual(h.wettedArea,34,'AbsTol',1e-9);
            t.verifyEqual(h.waterplaneArea,0);
            p.elevation=1.5; h=hydrostatics(t.Hull,p,t.Env);
            t.verifyEqual(h.volume,0);
            p.elevation=0; h=hydrostatics(t.Hull,p,t.Env);
            t.verifyEqual(h.volume,12,'AbsTol',1e-9);
        end
        function windingAndInvalidMeshes(t)
            hull=t.Hull; hull.faces=hull.faces(:,[1 3 2]);
            h=hydrostatics(hull,t.Pose,t.Env); t.verifyEqual(h.volume,4,'AbsTol',1e-9);
            hull=t.Hull; hull.faces(1,:)=hull.faces(1,[1 3 2]);
            t.verifyError(@()hydrostatics(hull,t.Pose,t.Env),'hull:Mesh');
            hull=t.Hull; hull.faces(1,:)=[];
            t.verifyError(@()hydrostatics(hull,t.Pose,t.Env),'hull:Mesh');
        end
        function draftAndImpossibleLoads(t)
            load=struct('mass',4000,'cg',[2 0 1]);
            r=solve_float(t.Hull,load,struct(),t.Env);
            t.verifyTrue(r.converged); t.verifyEqual(r.hydrostatics.draft,.5,'AbsTol',1e-9);
            t.verifyLessThan(abs(r.massResidual)/load.mass,1e-8);
            load.mass=12000;
            t.verifyError(@()solve_float(t.Hull,load,struct(),t.Env),'hull:Loading');
            load.mass=13000;
            t.verifyError(@()solve_float(t.Hull,load,struct(),t.Env),'hull:Loading');
            load.mass=-1;
            t.verifyError(@()solve_float(t.Hull,load,struct(),t.Env),'MATLAB:expectedPositive');
        end
        function analyticHeelAndSymmetry(t)
            a=10; p=t.Pose; p.heel=a; p.elevation=cosd(a);
            h=hydrostatics(t.Hull,p,t.Env);
            % Section wet depth .5 + y*tan(a), for -1 <= y <= 1.
            yB=(2/3*tand(a))/1;
            zB=1.25-tand(a)^2/3;
            t.verifyEqual(h.volume,4,'AbsTol',1e-9);
            t.verifyEqual(h.centreOfBuoyancyBody,[2 yB zB],'AbsTol',1e-9);
            t.verifyEqual(h.waterplaneArea,8/cosd(a),'RelTol',1e-9);
            t.verifyEqual(h.wettedArea,14,'AbsTol',1e-9);
            p.heel=-a; other=hydrostatics(t.Hull,p,t.Env);
            t.verifyEqual(other.volume,h.volume,'AbsTol',1e-9);
            t.verifyEqual(other.centreOfBuoyancyBody,h.centreOfBuoyancyBody.*[1 -1 1],'AbsTol',1e-9);
        end
        function rigidFrameEquivalence(t)
            p=t.Pose; p.heel=9; p.trim=-4;
            h=hydrostatics(t.Hull,p,t.Env);
            hull=t.Hull; hull.vertices=hull.vertices*hgeom.rotation(p).';
            p.heel=0; p.trim=0; other=hydrostatics(hull,p,t.Env);
            t.verifyEqual(other.volume,h.volume,'AbsTol',1e-9);
            t.verifyEqual(other.centreOfBuoyancyWorld,h.centreOfBuoyancyWorld,'AbsTol',1e-9);
            t.verifyEqual(other.wettedArea,h.wettedArea,'AbsTol',1e-9);
        end
        function analyticTrim(t)
            a=5; p=t.Pose; p.trim=a;
            % Put the waterline through body z=1 at midlength x=2.
            p.elevation=cosd(a)-2*sind(a);
            h=hydrostatics(t.Hull,p,t.Env);
            expectedX=2-(8/3)*tand(a);
            expectedZ=1.25-(4/3)*tand(a)^2;
            t.verifyEqual(h.volume,4,'AbsTol',1e-9);
            t.verifyEqual(h.centreOfBuoyancyBody,[expectedX 0 expectedZ],'AbsTol',1e-9);
            t.verifyEqual(h.waterplaneArea,8/cosd(a),'AbsTol',1e-9);
            t.verifyEqual(h.wettedArea,14,'AbsTol',1e-9);
        end
        function edgesAndVertices(t)
            p=t.Pose; p.heel=20;
            w=t.Hull.vertices*hgeom.rotation(p).'; p.elevation=w(4,3);
            h=hydrostatics(t.Hull,p,t.Env);
            t.verifyGreaterThan(h.volume,0); t.verifyLessThan(h.volume,12);
            p.trim=7; w=t.Hull.vertices*hgeom.rotation(p).'; p.elevation=w(4,3);
            h=hydrostatics(t.Hull,p,t.Env);
            t.verifyGreaterThan(h.volume,0); t.verifyLessThan(h.volume,12);
        end
        function stableAndUnstableEquilibrium(t)
            load=struct('mass',4000,'cg',[2 0 1]); opt=struct('mode','equilibrium');
            r=solve_float(t.Hull,load,opt,t.Env);
            t.verifyTrue(r.converged); t.verifyEqual(r.stability,'stable');
            opt.heel=3; opt.trim=2; r=solve_float(t.Hull,load,opt,t.Env);
            t.verifyTrue(r.converged);
            t.verifyEqual([r.pose.heel r.pose.trim],[0 0],'AbsTol',1e-5);
            t.verifyLessThan(norm(r.momentResidualWorld),.01);
            load.cg=[2 0 -1]; opt.heel=0; opt.trim=0;
            r=solve_float(t.Hull,load,opt,t.Env);
            t.verifyTrue(r.converged); t.verifyEqual(r.stability,'unstable');
            load.cg=[2 0 1.25-2/3]; r=solve_float(t.Hull,load,opt,t.Env);
            t.verifyEqual(r.stability,'near_neutral');
        end
        function reportsNonconvergence(t)
            load=struct('mass',4000,'cg',[2 .1 1]);
            opt=struct('mode','equilibrium','maxIterations',0);
            r=solve_float(t.Hull,load,opt,t.Env);
            t.verifyFalse(r.converged); t.verifyEqual(r.stability,'not_evaluated');
        end
        function offsetLoadAndBounds(t)
            load=struct('mass',4000,'cg',[2 .08 1]);
            opt=struct('mode','equilibrium');
            r=solve_float(t.Hull,load,opt,t.Env);
            t.verifyTrue(r.converged); t.verifyGreaterThan(r.pose.heel,0);
            t.verifyEqual(r.stability,'stable');
            t.verifyLessThan(norm(r.momentResidualWorld),.01);
            opt.angleBounds=[-1 1;-1 1];
            r=solve_float(t.Hull,load,opt,t.Env);
            t.verifyFalse(r.converged);
            t.verifyLessThanOrEqual(abs(r.pose.heel),1);
            opt.heel=2;
            t.verifyError(@()solve_float(t.Hull,load,opt,t.Env),'hull:Options');
        end
        function curvedRefinement(t)
            errors=[];
            for n=[16 32 64]
                a=(0:n-1)'*2*pi/n;
                s=struct('x',{0,4},'yz',{[cos(a) sin(a)],[cos(a) sin(a)]});
                [hull.vertices,hull.faces]=hgeom.section_mesh(s);
                p=t.Pose; p.elevation=0;
                h=hydrostatics(hull,p,t.Env);
                % Half-submerged circular cylinder, including semicircular end caps.
                errors(end+1)=abs(h.volume-2*pi)+abs(h.wettedArea-5*pi); %#ok<AGROW>
            end
            t.verifyLessThan(errors(2),errors(1)/3);
            t.verifyLessThan(errors(3),errors(2)/3);
        end
    end
end
