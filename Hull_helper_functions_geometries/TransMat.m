% TransMat
%
% Created by Anders Rosén, aro@kth.se.
%
% Function for defining a matrix for coordiate
% transformation corresponding to
% rotations (eta4,eta5,eta6) around the x-, y-, and z-axes
% in the following order:
% i)   eta6 [rad] rotation around the z-axis, "yaw";
% ii)  eta5 [rad] rotation around the y-axis, "pitch";
% iii) eta4 [rad] rotation around the x-axis, "roll".
%
% INPUT:
% eta6, eta5, eta4
%
% OUTPUT:
% T : transformation matrix

function [T] = TransMat(eta6,eta5,eta4)

% ----------------------------------
T(1,1) = cos(eta6)*cos(eta5);
T(1,2) = sin(eta6)*cos(eta5);
T(1,3) = -sin(eta5);
T(2,1) = cos(eta6)*sin(eta5)*sin(eta4)-sin(eta6)*cos(eta4);
T(2,2) = sin(eta6)*sin(eta5)*sin(eta4)+cos(eta6)*cos(eta4);
T(2,3) = cos(eta5)*sin(eta4);
T(3,1) = cos(eta6)*sin(eta5)*cos(eta4)+sin(eta6)*sin(eta4);
T(3,2) = sin(eta6)*sin(eta5)*cos(eta4)-cos(eta6)*sin(eta4);
T(3,3) = cos(eta5)*cos(eta4);
% ---------------------------------