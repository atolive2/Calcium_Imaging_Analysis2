function [ cat_bytrial, sum_byROI  ] = cat_response( vals_bytrial, threshold )
%cat_response categorizes ear trial for each ROI as response or no
%response.
%   inputs: vals_bytrial = numROIs by numtrials cell array of single values
%           threshold = scalar, minimum value that peak must pass in order
%           to be a response.
%   output: cat_bytrial = numROIs by numtrials logical array,
%           1 = val > threshold

data = cell2mat(vals_bytrial)

cat_bytrial = data > threshold
sum_byROI = sum(cat_bytrial,2)

end

