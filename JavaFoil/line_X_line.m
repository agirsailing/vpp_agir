function p = line_X_line(p1,p2,p3,p4);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function returns the intersection point p between the two 2D line
% segments p1-p2 and p3-p4. If the line segments do not intersect the
% function returns p=[]; Points are defined as vectors [x,y].
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x1  = p1(1); y1 = p1(2);
x2  = p2(1); y2 = p2(2);
x3  = p3(1); y3 = p3(2);
x4  = p4(1); y4 = p4(2);
p   = [];
n = (y4-y3)*(x2-x1)-(x4-x3)*(y2-y1);
if abs(n)>0.0; %(1e-25): 
	ua = ((x4-x3)*(y1-y3)-(y4-y3)*(x1-x3))/n;
	ub = ((x2-x1)*(y1-y3)-(y2-y1)*(x1-x3))/n;
    if ((ua>=0.0) & (ua<=1.0) & (ub>=0.0) & (ub<=1.0));
		% Intersection is within both line-segments
		x = x1+ua*(x2-x1);
		y = y1+ua*(y2-y1);
		p=[x,y];
    end
end