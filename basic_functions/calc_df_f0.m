function [ df_f0 ] = calc_df_f0( signal, stimulus_frame )
%calc_df_f0 takes signal and uses the pre-stimulus area to calculate df/f0
%for each frame over each ROI
%   inputs: signal = cell array with dims num_ROIs by num_trials, must be
%   final manipulated signal
%            stimulus_frame = frame number where stimulus starts
%   output: cell array with dims num_ROIs by num_trials

for i = 1:size(signal,1)
    for j = 1:size(signal,2)
        F0 = mean(signal{i,j}(1:stimulus_frame));
        df_f0{i,j} = (signal{i,j} - F0)/F0;
    end
end

end

