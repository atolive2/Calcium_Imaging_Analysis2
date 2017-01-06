function [ background ] = calc_background( data )
%calc_background takes tadpole.trial_split and generates a background value
%for each trial using last ROI.
%   input = data, a cell array of matrices of dims num_ROIs x num_trials
%   output = background, a vector of length num_trials. 

sizes = size(data)
for i = 1:size(data,2)
    background(i)=mean(data{end,i})
end

