function [ trial_split ] = split_into_trials( data, trial_length )
%split_into_trials take a type double matrix of F values for an entire experiment
%and split into individual trials. Requires all trials to be same length.
%   inputs: data = somaticF or neuropilF
%           trial_length = number of frames in a trial
%   output: trial_split = a cell array of dim num_trials x num_ROIs, each
%   array containing a vector of length trial_length.
if length(trial_length) == 1
    num_trials = size(data,2) / trial_length
    if floor(num_trials) ~= num_trials
        sprintf('error: num_trials is not an integer, value = %d', num_trials)
        return
    end
    frame_end = 0;
    frame_start = 1;
    split_trials=cell(size(data,1),num_trials);
    for i = 1:num_trials
        frame_end = frame_end + trial_length
        for j = 1:size(data,1)
        split_trials{j,i} = data(j, frame_start:frame_end);
        end
        frame_start = frame_end + 1
    end
    
else length(trial_length) ~= 1
    trials = [];
    for k = 1:12
        trials = [trials; trial_length]
    end
    all_trial_lengths = reshape(trials, 1,[])
    frame_end = 0;
    frame_start = 1;
    num_trials = length(all_trial_lengths)
    split_trials=cell(size(data,1),num_trials);
    for i = 1:num_trials
        frame_end = frame_end + all_trial_lengths(i)
        for j = 1:size(data,1)
        split_trials{j,i} = data(j, frame_start:frame_end);
        end
        frame_start = frame_end + 1
    end
end
trial_split = split_trials
end

% 
% frame_end = 0
% frame_start = 1
% frames = trial_length*num_trials
% for ii = 1:numtrialblocks;
%     frame_end = frame_end + frames
%     tadpole1.rawblock_somas{ii} = [raw_data(frame_start:frame_end,1:somaticRoiCounter)];
%     frame_start = frame_end + 1
% end