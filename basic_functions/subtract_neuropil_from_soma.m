function [ signal ] = subtract_neuropil_from_soma( soma_data, neuropil_data )
%subtract_neuropil_from_soma takes background subtracted somatic and
%neuropil F of the same size and subtracts to get signal
%   inputs: soma_data = cell array with dims num_ROIs by num_trials
%           neuropil_data = cell array with dims num_ROIs by num_trials
%   output: signal = cell array with dims num_ROIs by num_trials
for i = 1:size(soma_data,1)
    for j = 1:size(soma_data,2)
        signal{i,j} = soma_data{i,j} - neuropil_data{i,j} * 0.7;
    end
end
