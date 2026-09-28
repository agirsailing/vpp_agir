function hull = box_hull(L,B,H)
%BOX_HULL Independent closed box fixture, deck z=0, bottom z=H.
arguments
    L (1,1) double {mustBePositive,mustBeFinite} = 4
    B (1,1) double {mustBePositive,mustBeFinite} = 2
    H (1,1) double {mustBePositive,mustBeFinite} = 1.5
end
hull.name = 'Box'; hull.units = 'm';
hull.vertices = [0 -B/2 0;L -B/2 0;L B/2 0;0 B/2 0; ...
    0 -B/2 H;L -B/2 H;L B/2 H;0 B/2 H];
hull.faces = [1 3 2;1 4 3;5 6 7;5 7 8;1 2 6;1 6 5; ...
    2 3 7;2 7 6;3 4 8;3 8 7;4 1 5;4 5 8];
end
