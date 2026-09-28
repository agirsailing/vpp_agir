% Independent box validation and examples of all three solver modes.
hull=box_hull(4,2,1.5);
env=struct('water',struct('rho',1000),'g',9.81);
pose=struct('elevation',1,'heel',0,'trim',0);
hs=hydrostatics(hull,pose,env);
fprintf('Volume %.9g m^3; wetted area %.9g m^2\n',hs.volume,hs.wettedArea);
disp(hs.centreOfBuoyancyBody);
loading=struct('mass',4000,'cg',[2 0 1]);
upright=solve_float(hull,loading,struct('mode','draft'),env);
heeled=solve_float(hull,loading,struct('mode','draft','heel',10),env);
free=solve_float(hull,loading,struct('mode','equilibrium','heel',3,'trim',2),env);
disp(free.pose); disp(free.stability);
