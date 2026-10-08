function [] = plotFoils(mFoilCoords,fFoilCoords,alpha1,alpha2)

% Plot foil assembly
if isempty(mFoilCoords)==0
figure(2)
plot(mFoilCoords(:,1),mFoilCoords(:,2)); hold on
plot(fFoilCoords(:,1),fFoilCoords(:,2)); hold on
axis equal
grid on
plot([0 1],[0 0],'--r');
title(['Foil assembly, AoA Main = ',num2str(alpha1*180/pi),'°',', AoA Flap = ',num2str(alpha2*180/pi),'°'])
end