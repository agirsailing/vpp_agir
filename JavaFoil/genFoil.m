function [mFoilCoords fFoilCoords] = genFoil(main,flap,output,sclFlapVec,gap,mAngle,fAngle,mPivot,fPivot)

% Open files containing Foil coordinates
main = [pwd,'/Foils/',main,'.txt'];
flap  = [pwd,'/Foils/',flap,'.txt'];

fid1 = fopen(main);
fid2 = fopen(flap);
fid3 = fopen([pwd,'/JavaFoil/',output], 'w');

% Reads foil section coordinates from JavaFoil output and stores in matrix
% [x y]
A           = fscanf(fid1,'%c');
ind         = A==',';
A(ind)      = '.';
mFoilCoords = str2num(A);

A          = fscanf(fid2,'%c');
ind        = A==',';
A(ind)     = '.';
fFoilCoords = str2num(A);

% Scale flap
fFoilCoords = [fFoilCoords(:,1)*sclFlapVec+max(mFoilCoords(1,:))*(gap+1) fFoilCoords(:,2)*sclFlapVec];

% Scale to unit chord length
scale        = max(fFoilCoords(:,1));
fFoilCoords  = fFoilCoords./scale;
mFoilCoords  = mFoilCoords./scale;

% Define flap and main foil pivot
fPivot = fPivot*max(mFoilCoords(:,1));

% Rotate flap foil
T1 = [cos(fAngle) -sin(fAngle);
     sin(fAngle) cos(fAngle)];
fFoilCoords      = [(fFoilCoords(:,1)-fPivot) fFoilCoords(:,2)]*T1;
fFoilCoords(:,1) = fFoilCoords(:,1) + fPivot;

% Rotate complete foil
T2 = [cos(mAngle) -sin(mAngle);
     sin(mAngle) cos(mAngle)];
 
mFoilCoords = [(mFoilCoords(:,1)-mPivot) mFoilCoords(:,2)]*T2;
mFoilCoords(:,1) = mFoilCoords(:,1) + mPivot;
fFoilCoords = [(fFoilCoords(:,1)-mPivot) fFoilCoords(:,2)]*T2;
fFoilCoords(:,1) = fFoilCoords(:,1) + mPivot;
            
% Write coordinates to file read by Java script
fprintf(fid3, '%12.6f %12.6f\n', mFoilCoords');
fprintf(fid3, '%12.6f %12.6f\n', [9999.9 9999.9]);
fprintf(fid3, '%12.6f %12.6f\n', fFoilCoords');

% Close file fid3
fclose(fid1);
fclose(fid2);
fclose(fid3);
     
