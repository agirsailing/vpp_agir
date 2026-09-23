function R = rotation(pose)
% Active rotations: positive heel lowers starboard; positive trim raises bow.
validateattributes(pose.heel, {'numeric'}, {'scalar','real','finite'});
validateattributes(pose.trim, {'numeric'}, {'scalar','real','finite'});
c = cosd(pose.heel); s = sind(pose.heel);
C = cosd(pose.trim); S = sind(pose.trim);
R = [C 0 S;0 1 0;-S 0 C]*[1 0 0;0 c -s;0 s c];
end
