% TransShipfixedGlobal
%
% Created by Anders Rosén, aro@kth.se.
%
% Function for transformation of coordinates
% from ship fixed local coordinate system xyz
% to global coordinate system XYZ
% according to
% X = T'*x + r
% where
% r are translations (eta1,eta2,eta3) along the X-, Y-, and Z-axes
% T is the transformation matrix defined by TransMat.m
% for rotations (eta4,eta5,eta6) around the x-, y-, and z-axes
% in the following order
% i)   eta6 [rad] rotation around the z-axis, "yaw";
% ii)  eta5 [rad] rotation around the y-axis, "pitch";
% iii) eta4 [rad] rotation around the x-axis, "roll".
%
% INPUT:
% x,y,z : local ship fixed coordinates
% eta   : [eta1 eta2 eta3 eta4 eta5 eta6]
%
% OUTPUT:
% X,Y,Z : global coordinates

function [X,Y,Z] = TransShipfixedGlobal(x,y,z,eta)

r = eta(1:3); eta4 = eta(4); eta5 = eta(5); eta6 = eta(6);

[T] = TransMat(eta6,eta5,eta4);
T   = T';
X   = T(1,1)*x+T(1,2)*y+T(1,3)*z+r(1);
Y   = T(2,1)*x+T(2,2)*y+T(2,3)*z+r(2);
Z   = T(3,1)*x+T(3,2)*y+T(3,3)*z+r(3);


