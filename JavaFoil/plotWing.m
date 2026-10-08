function [] = plotWing(wing,mFoilCoords,fFoilCoords)

sections_xzc = wing.sections_xzc;

figure
for sec_=1:length(sections_xzc(:,1))-1
    x_flap = [];
    y_flap = [];
    z_flap = [];
    
    x_main = [];
    y_main = [];
    z_main = [];
    
    [~,yMainMin] = min(mFoilCoords(:,1));
    [~,yFlapMin] = min(fFoilCoords(:,1));
    for i_strip = 1:length(fFoilCoords(:,1))-1
        x_flap(:,i_strip) = [fFoilCoords(i_strip,1)*sections_xzc(sec_,3)+sections_xzc(sec_,1); fFoilCoords(i_strip+1,1)*sections_xzc(sec_,3)+sections_xzc(sec_,1); fFoilCoords(i_strip+1,1)*sections_xzc(sec_+1,3)+sections_xzc(sec_,1); fFoilCoords(i_strip,1)*sections_xzc(sec_+1,3)+sections_xzc(sec_,1);];
        y_flap(:,i_strip) = [fFoilCoords(i_strip,2)*sections_xzc(sec_,3)-mFoilCoords(yMainMin,2)*sections_xzc(sec_,3); fFoilCoords(i_strip+1,2)*sections_xzc(sec_,3)-mFoilCoords(yMainMin,2)*sections_xzc(sec_,3); fFoilCoords(i_strip+1,2)*sections_xzc(sec_+1,3)-mFoilCoords(yMainMin,2)*sections_xzc(sec_+1,3); fFoilCoords(i_strip,2)*sections_xzc(sec_+1,3)-mFoilCoords(yMainMin,2)*sections_xzc(sec_+1,3);];
        z_flap(:,i_strip) = [sections_xzc(sec_,2); sections_xzc(sec_,2); sections_xzc(sec_+1,2); sections_xzc(sec_+1,2)];
    end
    for i_strip = 1:length(mFoilCoords(:,1))-1
        x_main(:,i_strip) = [mFoilCoords(i_strip,1)*sections_xzc(sec_,3)+sections_xzc(sec_,1); mFoilCoords(i_strip+1,1)*sections_xzc(sec_,3)+sections_xzc(sec_,1); mFoilCoords(i_strip+1,1)*sections_xzc(sec_+1,3)+sections_xzc(sec_,1); mFoilCoords(i_strip,1)*sections_xzc(sec_+1,3)+sections_xzc(sec_,1);];
        y_main(:,i_strip) = [mFoilCoords(i_strip,2)*sections_xzc(sec_,3)-mFoilCoords(yMainMin,2)*sections_xzc(sec_,3); mFoilCoords(i_strip+1,2)*sections_xzc(sec_,3)-mFoilCoords(yMainMin,2)*sections_xzc(sec_,3); mFoilCoords(i_strip+1,2)*sections_xzc(sec_+1,3)-mFoilCoords(yMainMin,2)*sections_xzc(sec_+1,3); mFoilCoords(i_strip,2)*sections_xzc(sec_+1,3)-mFoilCoords(yMainMin,2)*sections_xzc(sec_+1,3)];
        z_main(:,i_strip) = [sections_xzc(sec_,2); sections_xzc(sec_,2); sections_xzc(sec_+1,2); sections_xzc(sec_+1,2)];
    end
patch(x_flap,y_flap,z_flap,ones(size(z_flap)),'FaceColor',[0 0 1],'FaceAlpha',0.3,'LineStyle','none'); hold on
plot3(mFoilCoords(:,1)*sections_xzc(sec_,3)+sections_xzc(sec_,1),mFoilCoords(:,2)*sections_xzc(sec_,3)-mFoilCoords(yMainMin,2)*sections_xzc(sec_,3),ones(size(mFoilCoords(:,1)))*sections_xzc(sec_,2),'k','linewidth',1.5)
patch(x_main,y_main,z_main,ones(size(z_main)),'FaceColor',[0 1 0],'FaceAlpha',0.3,'LineStyle','none')
plot3(fFoilCoords(:,1)*sections_xzc(sec_,3)+sections_xzc(sec_,1),fFoilCoords(:,2)*sections_xzc(sec_,3)-mFoilCoords(yMainMin,2)*sections_xzc(sec_,3),ones(size(fFoilCoords(:,1)))*sections_xzc(sec_,2),'k','linewidth',1.5)
axis equal

end
plot3(mFoilCoords(:,1)*sections_xzc(sec_+1,3)+sections_xzc(sec_+1,1),mFoilCoords(:,2)*sections_xzc(sec_+1,3)-mFoilCoords(yMainMin,2)*sections_xzc(sec_+1,3),ones(size(mFoilCoords(:,1)))*sections_xzc(sec_+1,2),'k','linewidth',1.5)
plot3(fFoilCoords(:,1)*sections_xzc(sec_+1,3)+sections_xzc(sec_+1,1),fFoilCoords(:,2)*sections_xzc(sec_+1,3)-mFoilCoords(yMainMin,2)*sections_xzc(sec_+1,3),ones(size(fFoilCoords(:,1)))*sections_xzc(sec_+1,2),'k','linewidth',1.5)
axis([-1 10 -3 3 0 38])
grid on
box on
view([10 38])
xlabel('X [m]')
ylabel('Y [m]')
zlabel('Z [m]')
