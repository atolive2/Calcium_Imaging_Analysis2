function [ area_bytrial ] = calc_area( signal, stimulus_start_frame )
%calc_area uses trapz to generate the area under the curve for a trial.
%   inputs: data = df/f0, a cell array with dims num_ROIs by num_trials
%           stimulus_start_frame = frame number where response starts
%   output: area_bytrial is a type double with dims num_ROIs by num_trials,
%   1 value per trial per ROI

for i = 1:size(signal,1)
    for j = 1:size(signal,2)
        area_bytrial{i,j}=trapz(signal{i,j}(stimulus_start_frame:end));
    end
end
end

