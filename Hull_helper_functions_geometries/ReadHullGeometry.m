% ReadHullGeometry
%
% Created by Anders Rosén, aro@kth.se.
%
% Function for reading a hull geometry defined in a britfair ".bri" text file
% into Matlab's work space.
%
% Called by either "[ship]=ReadHullGeometry"
% where you will be asked to choose a hull geometry file,
% or by "[ship]=ReadHullGeometry(fname)"
% where you should specify the hull geometry file name in fname,
% for example fname = 'Box.bri' or fname = 'AxelJohnson_Container.bri'
%
% The output is the structured array "ship" which includes:
%   ship.x; ship.y; ship.z : hull offset point coordinates in a ship fixed coord syst
%   ship.ns : number of sections
%   ship.np : number of offset points for each section

%%
function [ship]=ReadHullGeometry(fname)

if nargin==0
    [fname,pathname]=uigetfile('*.bri','Choose hull geometry file!');
else
    pathname = 'HullGeometries\'; %mac-users, replace '\' by '/'!
end

[x0,y0,z0,np0,K,namn]=ReadBritfair2(pathname,fname);
ns = size(x0,2); ns0 = ns;
np = 2*np0+1;
x=zeros(max(np),ns0); y=zeros(max(np),ns0); z=zeros(max(np),ns0);
for I=1:ns
    x(1:np(I),I) = [x0(1:np0(I),I) ; x0(np0(I),I) ; flipud(x0(1:np0(I),I))]';
    y(1:np(I),I) = [y0(1:np0(I),I) ; 0 ; flipud(-y0(1:np0(I),I))]';
    z(1:np(I),I) = [z0(1:np0(I),I) ; z0(np0(I),I) ; flipud(z0(1:np0(I),I))]';
end

ship.geometryfile=fname; ship.x=x; ship.y=y; ship.z=z; ship.np=np; ship.ns=ns;
state.geometryfile = fname;

return;


% ----- ReadBritfair2 ------------------------------
function [x,y,z,np,K,namn]=ReadBritfair2(pathname,fname)
fid=fopen([pathname fname],'r');
namn = fgetl(fid); slaskrad = fgetl(fid);
[SecHead,count]=fscanf(fid,'%f',3);
J = 0;
while SecHead(1)>0
    if SecHead(3)>-999
        J      = J+1; I = 1;
        K(I,J) = 0;
        np(J)  = SecHead(1);
        Q      = 1:np(J);
        xk     = SecHead(3);
    else
        K(I,J) = 1;
        np(J)  = np(J)+SecHead(1)-1;
        Q      = I:(I+SecHead(1)-1);
        xk     = x(1,J);
    end
    for I = Q
        [yz,count]=fscanf(fid,'%f',2) ;
        x(I,J) = xk;  y(I,J) = yz(1);  z(I,J) = yz(2);
    end
    [Slask,count]=fscanf(fid,'%f',1) ;
    if Slask~=0
        [Slask,count]=fscanf(fid,'%f',2) ;
    end
    [SecHead,count]=fscanf(fid,'%f',3) ;
end
fclose(fid);
return;
% -----------------------------------------------------
