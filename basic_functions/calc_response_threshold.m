function [ cutoff, std_dev, avg_alltrial ] = calc_response_threshold( data, desired_cutoff_sd )
%calc_response_threshold determines the threshold for response. 
%   inputs: data = tadpole{i}.datavalues (e.g. meanpeak_bytrial) 
%           desired cutoff = 1, 2 or 3 SD above mean
%   output: cutoff, a scalar value representing the threshold.

%if class(data) == cell
    data = cell2mat(data)
%end

datav = reshape(data,[],1);
avg_alltrial = mean(datav);
std_dev = std(datav);
cutoff = avg_alltrial + (std_dev * desired_cutoff_sd);

end

