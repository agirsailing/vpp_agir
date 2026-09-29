function [foil_data] = extract_data(table,Re)
% takes the xfoil data and keeps only the values corresponding to the
% reynoldsnumber


% sort on Reynoldsnumber in ascending order
T = sortrows(table,"Re","ascend");

% find ALL indicies for the closest Re values in the table T
[~, idx] = min(abs(T.Re - Re));

% Select all rows corresponding to the closest Reynolds number
closestRe = T.Re(idx);

% Extract all rows matching the closest Reynolds number
foil_data = T(T.Re == closestRe, :);

end