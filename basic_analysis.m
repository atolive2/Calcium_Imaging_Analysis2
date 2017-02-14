%% basic properties (for committee meeting 2/22/2017)

% Percent response by cell (all stims combined)
t=1
all_prop_responses = [];
for t = 1:length(tadpole)
    stimuli = tadpole{1,t}.stimorder ~= 4
    total_stimuli = sum(stimuli)
    prop_responses = zeros(size(tadpole{1,t}.boolean_response,1),1);
    for i=1:length(stimuli)
        if stimuli(i)
            for j = 1:size(tadpole{1,t}.boolean_response,1)
                prop_responses(j,1) = prop_responses(j,1) + tadpole{1,t}.boolean_response(j,i);
            end
        end
    end
    prop_responses(:,2) = prop_responses(:,1)/total_stimuli;
    all_prop_responses = [all_prop_responses; prop_responses];
end
n = size(all_prop_responses,1)
hist(all_prop_responses(:,2),25)
title('Proportion of stimuli responded to by cell, all cells, all stims')
figure;
hist(all_prop_responses(:,1),25)
title('Number of stimuli responded to by cell, all cells, all stims')

% percent response by cell by condition
