function [ background_sub ] = subtract_background( data, background_vals )
%subtract_background takes trial split data and subtracts the background
%for that trial from every frame for every ROI
%   inputs: data = tadpole.trial_split
%           background_vals = tadpole.background
%   output: background_sub = cell array with dims num_ROIs by num_trials
%   that is raw - background

for i = 1:size(data,1)
    for j = 1:size(data,2)
        background_sub{i,j} = data{i,j} - background_vals(j);
    end
end

end

