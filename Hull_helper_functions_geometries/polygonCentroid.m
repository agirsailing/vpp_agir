%POLYGONCENTROID compute centroid (center of mass) of a polygon
%
%   PT = polygonCentroid(POINTS)
%   PT = polygonCentroid(PTX, PTY)
%   Computes center of mass of a polygon defined by POINTS. POINTS is a
%   [N*2] array of double.
function [pt As] = polygonCentroid(varargin)
%   See also:
%   polygons2d, polygonArea, drawPolygon
%
%   ---------
%   author : David Legland 
%   INRA - TPV URPOI - BIA IMASTE
%   created the 05/05/2004.

if nargin==1
    var = varargin{1};
    px = var(:,1);
    py = var(:,2);
elseif nargin==2
    px = varargin{1};
    py = varargin{2};
end

% Algorithme P. Bourke
sx = 0;
sy = 0;
N = length(px);
for i=1:N-1
    sx = sx + (px(i)+px(i+1))*(px(i)*py(i+1) - px(i+1)*py(i));
    sy = sy + (py(i)+py(i+1))*(px(i)*py(i+1) - px(i+1)*py(i));
end
sx = sx + (px(N)+px(1))*(px(N)*py(1) - px(1)*py(N));
sy = sy + (py(N)+py(1))*(px(N)*py(1) - px(1)*py(N));
As = polygonArea(px, py);
pt = [sx sy]/6/As;