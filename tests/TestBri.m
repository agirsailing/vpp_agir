classdef TestBri < matlab.unittest.TestCase
    properties
        File
        Options
    end
    methods(TestClassSetup)
        function projectPath(t)
            t.applyFixture(matlab.unittest.fixtures.PathFixture(fileparts(fileparts(mfilename('fullpath')))));
        end
    end
    methods(TestMethodSetup)
        function setup(t)
            t.File=[tempname '.bri'];
            fid=fopen(t.File,'w'); fclose(fid);
            t.addTeardown(@()delete(t.File));
            t.Options=struct('lengthScale',1,'xHeaderColumn',2, ...
                'axisSigns',[1 1 1],'originOffset',[0 0 0],'halfHull',false);
        end
    end
    methods(Test)
        function fullBoxWithDuplicates(t)
            t.writeBox(false,true);
            hull=read_bri(t.File,t.Options);
            t.verifyEqual(hull.source.header,17);
            t.verifyEqual(hull.sections(1).sourceTrailer,42);
            t.verifyEqual(size(hull.sections(1).yz,1),4);
            env=struct('water',struct('rho',1000),'g',9.81);
            hs=hydrostatics(hull,struct('elevation',1,'heel',0,'trim',0),env);
            t.verifyEqual(hs.volume,4,'AbsTol',1e-9);
            t.verifyEqual(hs.wettedArea,14,'AbsTol',1e-9);
        end
        function mirroredBox(t)
            t.writeBox(true,false); opts=t.Options; opts.halfHull=true;
            hull=read_bri(t.File,opts);
            env=struct('water',struct('rho',1000),'g',9.81);
            hs=hydrostatics(hull,struct('elevation',1,'heel',0,'trim',0),env);
            t.verifyEqual(hs.volume,4,'AbsTol',1e-9);
            t.verifyEqual(hs.wettedArea,14,'AbsTol',1e-9);
        end
        function explicitMapping(t)
            t.writeBox(false,false); o=t.Options;
            o.lengthScale=2; o.axisSigns=[-1 -1 -1]; o.originOffset=[8 0 3];
            hull=read_bri(t.File,o);
            t.verifyEqual([hull.sections.x],[0 8]);
            t.verifyEqual(min(hull.vertices),[0 -2 0],'AbsTol',1e-9);
            t.verifyEqual(max(hull.vertices),[8 2 3],'AbsTol',1e-9);
            o.xHeaderColumn=3; o.axisSigns=[1 1 1]; o.originOffset=[0 0 0];
            hull=read_bri(t.File,o); t.verifyEqual([hull.sections.x],[20 28]);
            o=rmfield(o,'halfHull');
            t.verifyError(@()read_bri(t.File,o),'hull:Options');
        end
        function malformedInput(t)
            fid=fopen(t.File,'w'); fprintf(fid,'Box\n1\n4 0 0\n0 0\n'); fclose(fid);
            t.verifyError(@()read_bri(t.File,t.Options),'hull:BRI');
            fid=fopen(t.File,'w'); fprintf(fid,'Box\n1\n4 0 0\nNaN 0\n0 1\n1 1\n1 0\n0\n0 0 0\n'); fclose(fid);
            t.verifyError(@()read_bri(t.File,t.Options),'hull:BRI');
            t.writeBox(false,false);
            text=fileread(t.File); text=strrep(text,'4 4 14','4 0 14');
            fid=fopen(t.File,'w'); fprintf(fid,'%s',text); fclose(fid);
            t.verifyError(@()read_bri(t.File,t.Options),'hull:BRI');
        end
        function selfIntersection(t)
            p=[-1 0;1 1.5;-1 1.5;1 0];
            s=struct('x',{0,4},'yz',{p,p});
            t.verifyError(@()hgeom.section_mesh(s),'hull:Sections');
        end
        function varyingCountsRetainCorners(t)
            p=[-1 0;1 0;1 1.5;-1 1.5]; q=[-1 0;0 0;1 0;1 1.5;-1 1.5];
            s=struct('x',{0,4},'yz',{p,q});
            [v,~,s]=hgeom.section_mesh(s);
            for k=1:numel(s)
                for j=1:size(s(k).yz,1)
                    point=[s(k).x s(k).yz(j,:)];
                    t.verifyLessThan(min(vecnorm(v-point,2,2)),1e-10);
                end
            end
        end
        function explicitSealing(t)
            fid=fopen(t.File,'w'); fprintf(fid,'Open\n1\n');
            for x=[0 4]
                fprintf(fid,'4 %g %g\n0 0\n0 0.1\n1 0.1\n1 1.5\n0\n',x,x);
            end
            fprintf(fid,'0 0 0\n'); fclose(fid);
            o=t.Options; o.halfHull=true;
            t.verifyError(@()read_bri(t.File,o),'hull:BRI');
            o.halfHullClosure='seal'; h=read_bri(t.File,o);
            t.verifyTrue(all([h.sections.closureAdjusted]));
            t.verifyEqual(h.sections(1).sourcePoints,[0 0;0 .1;1 .1;1 1.5]);
        end
        function provisionalKayakRejectsTouchingContour(t)
            root=fileparts(fileparts(mfilename('fullpath')));
            o=t.Options; o.halfHull=true; o.axisSigns=[1 1 -1];
            o.halfHullClosure='seal';
            % Station 1.4 returns to y=0 internally, then extends out again.
            % Mirroring creates intersecting/touching contours; do not repair it.
            t.verifyError(@()read_bri(fullfile(root,'data','kayak.bri'),o),'hull:Sections');
        end
    end
    methods(Access=private)
        function writeBox(t,half,duplicate)
            if half, p=[0 0;1 0;1 1.5;0 1.5];
            else, p=[-1 0;1 0;1 1.5;-1 1.5]; end
            if duplicate, p=[p(1,:);p;p(1,:)]; end
            fid=fopen(t.File,'w'); cleanup=onCleanup(@()fclose(fid)); %#ok<NASGU>
            fprintf(fid,'Box\n17\n');
            for x=[0 4]
                fprintf(fid,'%d %g %g\n',size(p,1),x,x+10);
                fprintf(fid,'%g %g\n',p.'); fprintf(fid,'42\n');
            end
            fprintf(fid,'0 0 0\n');
        end
    end
end
