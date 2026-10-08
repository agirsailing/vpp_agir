function [cl cd mFoilCoords fFoilCoords] = calc_clcd(wing,alpha1,alpha2,twAngle,flag)

% Foil configuration parameters, do not change.
gap     = wing.gap;
fPivot  = wing.fPivot;
mPivot  = wing.mPivot;
main    = wing.main;           % Main foil
flap    = wing.flap;           % Flap foil
sclFlap = wing.sclFlap;        % Relation between c1 and c2

output = 'Wing1.txt';                   % Resulting 2-foil section
Re     = 8e6;                           % Reynolds number, assume constant.

if flag==1
    for i = 1:length(alpha1)
        for ii = 1:length(alpha2)

            % Generate foil combination
            [mFoilCoords fFoilCoords] = genFoil(main,flap,output,sclFlap,gap,alpha1(i),alpha2(ii),mPivot,fPivot);

            % Calculate lift and drag coeffs
            [cl(i,ii), cd(i,ii)] = FoilCalc(output,Re);
        end
    end
    disp('precalculations done')
    % Smooth cd for equation solving
    sz = size(cd);
    if mod(sz(1),2)==0; sz(1)=sz(1)-1; end
    if mod(sz(2),2)==0; sz(2)=sz(2)-1; end
    cd_     = smoothn(cd,sz,'gaussian',[0.15 0.15]);
    cl_     = cl;
    alpha1_ = alpha1;
    alpha2_ = alpha2;
    save('preCalc.mat','alpha1_','alpha2_','cl_','cd_');
    return
else  
    % Here interpolation takes place if flag=0;
    load(wing.preCalc)
    if length(alpha1)>1
        % Interpolate for cl and cd matrices at twist angle twAngle
        [alpha1_int alpha2_int] = meshgrid(alpha1_,alpha2_);
        cl = interp2(alpha1_int,alpha2_int,cl_',alpha1',alpha2,'cubic')';
        cd = interp2(alpha1_int,alpha2_int,cd_',alpha1',alpha2,'cubic')';
    else
        [alpha1_int alpha2_int] = meshgrid(alpha1_,alpha2_);
        cl = interp2(alpha1_int,alpha2_int,cl_',alpha1'-twAngle/8,alpha2-twAngle,'cubic')';
        cd = interp2(alpha1_int,alpha2_int,cd_',alpha1'-twAngle/8,alpha2-twAngle,'cubic')';
    end
    mFoilCoords = [];
    fFoilCoords = [];
end
