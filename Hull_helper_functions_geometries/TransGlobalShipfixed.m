% TransGlobalShipfixed
%
% Created by Anders Rosén, aro@kth.se.
%
% Function for transformation of coordinates
% from global coordinate system XYZ
% to ship fixed local coordinate system xyz
% according to
% x=T*(X-r)
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
% X,Y,Z : global coordinates
% eta   : [eta1 eta2 eta3 eta4 eta5 eta6]
%
% OUTPUT:
% x,y,z : local ship fixed coordinates

function [x,y,z] = TransGlobalShipfixed(X,Y,Z,eta)

r = eta(1:3); eta4 = eta(4); eta5 = eta(5); eta6 = eta(6);

[T] = TransMat(eta6,eta5,eta4);
x   = T(1,1)*(X-r(1))+T(1,2)*(Y-r(2))+T(1,3)*(Z-r(3));
y   = T(2,1)*(X-r(1))+T(2,2)*(Y-r(2))+T(2,3)*(Z-r(3));
z   = T(3,1)*(X-r(1))+T(3,2)*(Y-r(2))+T(3,3)*(Z-r(3));


