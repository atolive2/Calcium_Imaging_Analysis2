function [ sorted ] = sort_ROIs( order, data )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% order = 1xlength(respROI)
% data = length(respROI) x length(respROI)
% sorted = data, but rows and cols are in the order of order

for dim1 = 1:length(order)
    for dim2 = 1:length(order)
        sorted(dim1, dim2) = data(order(dim1), order(dim2));
    end
end

end

