function [ stddev ] = std_by_stimtype ( bytrial, stimmask )
%mean_by_stimtype takes the by trial data and the stim order and returns average values for each ROI. 
%   input: bytrial= type cell array with dims num_ROIs by num_trials
%          stimorder = type double, identifies stim type
%          boolean_response = include in the avg
%   output: avg = mean of all presentations of a given stimulus.

% this is a copy of mean_by_stimtype2 that does Standard deviation instead

for i = 1:size(stimmask,2) % over each stim type
    for j = 1:size(bytrial,1) %over each ROI

        
        stddev(i,j) = std(bytrial(j, stimmask(:,i)));
    end
end
   
end

