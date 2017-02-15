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

for t = 1:length(tadpole)
    conds= unique(tadpole{1,t}.stimorder);
    for i = 1:length(conds)
        trial_count(conds(i),t) = sum(tadpole{1,t}.stimorder(:) == conds(i));
    end
    count_resp_bycond = zeros(size(tadpole{1,t}.boolean_response, 1), length(conds));
    for i=1:length(conds) %over all stimulus conditions
        for j = 1:size(tadpole{1,t}.stimmask,1) %over all trials
            if tadpole{1,t}.stimmask(j,i) %1=this trial presented stimulus i
                for k = 1:size(tadpole{1,t}.boolean_response, 1)
                    count_resp_bycond(k,i) = count_resp_bycond(k,i) + tadpole{1,t}.boolean_response(k,j);
                end
            end
        end
    end
    count_resp_bycond_all{1,t} = count_resp_bycond;
    clear('conds')
end

% generate histograms by stim type
histdata=[];
for t = 1:length(tadpole)
    histdata=[histdata; count_resp_bycond_all{1,t}(:,1:4)];
end

for i = 1:size(histdata,2)
    figure;
    hist(histdata(:,i),25)
    title(sprintf('All response counts all presentations of stim %d', i))
end

