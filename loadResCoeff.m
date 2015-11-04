function [coeff_res, iterTime, objectValue] = loadResCoeff(fileName, headlines)
coeff_res = [];
iterTime = [];
objectValue = [];

fid = fopen(fileName, 'r');
cells = textscan(fid, '%f %f %f %f %f %f %f %f %f %f %f %f','HeaderLines',headlines, 'EmptyValue',0);
coeff_res = [cells{1}(2:2:end),cells{2}(2:2:end),cells{3}(2:2:end),cells{4}(2:2:end),cells{5}(2:2:end),cells{6}(2:2:end),cells{7}(2:2:end),cells{8}(2:2:end),cells{9}(2:2:end),cells{10}(2:2:end),cells{11}(2:2:end),cells{12}(2:2:end)];
iterTime = cells{1}(1:2:end);
objectValue = cells{2}(1:2:end);
% coeff_res = arrayfun(@(x));

fclose(fid);
end