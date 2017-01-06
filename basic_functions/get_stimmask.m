function [ stimmask ] = get_stimmask( stimorder )
%get_stimmask takes stim order and makes a column for each value of stimorder with 1's where true. 
%   input: stimorder = vector of length num_trials
%   output: stimmask = logical of size num_trials by unique(stimorder). 
%           1 = that stim was presented
%           0 = that stim was not presented

C = unique(stimorder)
stimmask = zeros(length(stimorder), length(C));
for i = 1:length(C)
    idxs = find(stimorder == C(i));
    single_stimmask=zeros(length(stimorder),1);
    for j = 1:length(idxs)
        single_stimmask(idxs) = 1;
    end
     stimmask(:,i) = single_stimmask;
     clear('idxs')
end

stimmask=logical(stimmask);
end
