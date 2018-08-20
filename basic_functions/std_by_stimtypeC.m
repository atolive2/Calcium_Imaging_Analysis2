function [ std_dev ] = std_by_stimtypeC ( data, include )
%mean_by_stimtype takes the by trial data and the stim order and returns average values for each ROI. 
%   input: bytrial= type double array with dims num_ROIs by num_trials
%          data = type double, identifies stim type
%          include = include in the avg, size num_stims x num_rois x
%          num_trials
%   output: avg = mean of all "good" presentations of a given stimulus, size num_stimtypes x num_rois.

for i = 1:size(include,1) % over each stim type
    for j = 1:size(data,1) %over each ROI
        
        
        std_dev(i,j) = std(data(j, include(i,j,:)));
    end
end
   
end