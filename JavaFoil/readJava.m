function [cl cd] = readJava(filepath)

% Open results file from JavaFoil
fid = fopen(filepath,'r');
for i=1:5
        fgets(fid);
end

% Read entire file
i = 0;
marker = fgets(fid);
while length(marker) > 1
    i = i + 1;
    if isempty(sscanf(marker,'%f')'); break; end
    out(i,:) = sscanf(marker,'%f')';
    marker = fgets(fid);
end
fclose(fid);

% Extract cl and cd and return to main
cl = out(:,2)';
cd = out(:,3)';

