function [Ywet,Zwet] = cut_wet_sections(Y,Z,H)
% Take out the wet secions


wet_points = Z<H; % All points (0=in air, 1=in water)

n = [];
m = [];

for k=2:length(wet_points)
    if wet_points(k) == wet_points(k-1)
        continue;
    else
        new_point = line_X_line([Y(k),Z(k)],[Y(k-1),Z(k-1)],[-100,H],[100,H]);
            if wet_points(k-1)==0 && wet_points(k)==1
                Y(k-1) = new_point(1);
                Z(k-1) = new_point(2);
                n = k-1;
            else if wet_points(k)==0 && wet_points(k-1)==1
                Y(k) = new_point(1);
                Z(k) = new_point(2);
                m = k;
                end
            end
    end
end

if isempty(n) % If the foil never enters the water
    n = 1;
end
if isempty(m) % If the foil never leaves the water
    m = length(Y);
end

Ywet = Y(n:m);
Zwet = Z(n:m);

% plot(Y,Z,'-ok','LineWidth',2), hold on
% % Add waterline to plot
% Zw = H*ones(1,length(Y));
% plot(Y,Zw,'LineWidth',2)
% plot(Ywet,Zwet,'-rs','LineWidth',2)
% axis equal

end

