%% basic properties (for committee meeting 2/22/2017)

%% Proportion responding (using boolean_response)

% Percent response by cell (all stims combined)
t=1
all_prop_responses = [];
for t = 1:length(tadpole)
    stimuli = tadpole{1,t}.stimorder ~= 4
    total_stimuli = sum(stimuli)
    prop_responses = zeros(size(tadpole{1,t}.boolean_response_sm,1),1);
    for i=1:length(stimuli)
        if stimuli(i)
            for j = 1:size(tadpole{1,t}.boolean_response_sm,1)
                prop_responses(j,1) = prop_responses(j,1) + tadpole{1,t}.boolean_response_sm(j,i);
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
    count_resp_bycond = zeros(size(tadpole{1,t}.boolean_response_sm, 1), length(conds));
    for i=1:length(conds) %over all stimulus conditions
        for j = 1:size(tadpole{1,t}.stimmask,1) %over all trials
            if tadpole{1,t}.stimmask(j,i) %1=this trial presented stimulus i
                for k = 1:size(tadpole{1,t}.boolean_response_sm, 1)
                    count_resp_bycond(k,i) = count_resp_bycond(k,i) + tadpole{1,t}.boolean_response_sm(k,j);
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

% Calculate proportion with response
for t = 1:length(tadpole)
    conds= unique(tadpole{1,t}.stimorder);
    for i=1:length(conds)
        prop_resp_bycond{1,t}(:,i) = count_resp_bycond_all{1,t}(:,i) / trial_count(conds(i),t);
    end
end

% hist proportion with response
prop_histdata=[];
for t = 1:length(tadpole)
    prop_histdata=[prop_histdata; prop_resp_bycond{1,t}(:,1:4)];
end

for i = 1:size(prop_histdata,2)
    figure;
    hist(prop_histdata(:,i),25)
    title(sprintf('All proportion responses all presentations of stim %d', i))
    axis([0 1 -inf inf])
end

%% Peak df/f0 values

% Gather all meanpeak values (peak +/-1 averaged) for all ROIs, all trials
peak_allhistdata = [];
for t = 1:length(tadpole)
    num_data = size(tadpole{1,t}.meanpeak_bytrial_sm, 1) * size(tadpole{1,t}.meanpeak_bytrial_sm, 2);
    trial_data = reshape(cell2mat(tadpole{1,t}.meanpeak_bytrial_sm), 1, num_data);
    peak_allhistdata = [peak_allhistdata trial_data];
    clear('trial_data', 'num_data')
end

outliers = find(peak_allhistdata > 30)
peak_allhistdata_nooutliers = peak_allhistdata(peak_allhistdata<30);
figure;
hist(peak_allhistdata)
hist(peak_allhistdata_nooutliers)
peak_allhistdata_small = peak_allhistdata(peak_allhistdata < 1);
peak_allhistdata_small = peak_allhistdata_small(peak_allhistdata_small>-0.2)
hist(peak_allhistdata_small,100)

%% How many cells per tadpole respond?
%t = 1

for t = 1:length(tadpole)
    total_trials = size(tadpole{1,t}.boolean_response_sm,2)
    total_resp = sum(tadpole{1,t}.boolean_response_sm,2)
    responders50 = total_resp > (total_trials / 2) 
    prop_responders50(t) = sum(responders50) / size(tadpole{1,t}.boolean_response_sm,1)
    responders25 = total_resp > (total_trials / 4) 
    prop_responders25(t) = sum(responders25) / size(tadpole{1,t}.boolean_response_sm,1)
    clear('responders50', 'responders25', 'total_trials', 'total_resp') 
end

bar(prop_responders25)
figure;
bar(prop_responders50)

% What is the primary modality of responders?
% pie chart
%t = 1
for t = 1:length(tadpole)
    for k = 1:size(tadpole{1,t}.peak_avg_sm, 2) % over all rois
        [tadpole{1,t}.unimax_peakavg_sm(1,k) tadpole{1,t}.unimax_stimtype_sm(1,k)] = max(tadpole{1,t}.peak_avg_sm(2:3, k));
        tadpole{1,t}.multimax_peakavg_sm(1,k) = tadpole{1,t}.peak_avg_sm(1, k);
    end
end


for t = 1:length(tadpole)
    total_trials = size(tadpole{1,t}.boolean_response_sm,2);
    total_resp = sum(tadpole{1,t}.boolean_response_sm,2);
    responders50 = total_resp > (total_trials / 2) ;
    responders25 = total_resp > (total_trials / 4) 
    primary_modality = [];
    for i = 1:length(responders50)
        if responders50(i,1)
            test = tadpole{1,t}.unimax_peakavg_sm(1,i) > tadpole{1,t}.multimax_peakavg_sm(1,i)
            if test
                primary_modality = [primary_modality; (tadpole{1,t}.unimax_stimtype_sm(1,i) +1)]
            else
                primary_modality = [primary_modality; 1]
            end
        end
    end
    primary_modality50{1,t} = primary_modality
    
    primary_modality = [];
    for i = 1:length(responders25)
        if responders25(i,1)
            test = tadpole{1,t}.unimax_peakavg_sm(1,i) > tadpole{1,t}.multimax_peakavg_sm(1,i)
            if test
                primary_modality = [primary_modality; (tadpole{1,t}.unimax_stimtype_sm(1,i) +1)]
            else
                primary_modality = [primary_modality; 1]
            end
        end
    end
    primary_modality25{1,t} = primary_modality
    %clear('primary_modality')
    clear('responders50', 'responders25', 'total_trials', 'total_resp')
end

for t = 1:length(tadpole)
    piect(1,1) = sum(primary_modality25{1,t} == 1)
    piect(1,2) = sum(primary_modality25{1,t} == 2)
    piect(1,3) = sum(primary_modality25{1,t} == 3)
    figure;
    labels = {'multi', 'vis', 'mech'}
    pie(piect, labels)
    title(sprintf('tadpole %d', t))
    fig_filename=sprintf('piechart_tadpole%d', t);
    saveas(gcf,fig_filename,'png');
    close;
    clear('fig_filename')
end

for t = 1:length(tadpole)
    piect(t,1) = sum(primary_modality25{1,t} == 1);
    piect(t,2) = sum(primary_modality25{1,t} == 2);
    piect(t,3) = sum(primary_modality25{1,t} == 3);
end
grandpie25 = sum(piect)
clear('piect')
for t = 1:length(tadpole)
    piect(t,1) = sum(primary_modality50{1,t} == 1)
    piect(t,2) = sum(primary_modality50{1,t} == 2)
    piect(t,3) = sum(primary_modality50{1,t} == 3)
end
grandpie50 = sum(piect)

figure;
labels = {'multi', 'vis', 'mech'}
pie(grandpie25, labels)
title('25% responders')
fig_filename='grandpie25_sm';
saveas(gcf,fig_filename,'png');

figure;
labels = {'multi', 'vis', 'mech'}
pie(grandpie50, labels)
title('50% responders')
fig_filename='grandpie50_sm';
saveas(gcf,fig_filename,'png');

%% Scatterplot of avg peak multi vs uni

for t = 1:length(tadpole)
    scatter(tadpole{1,t}.unimax_peakavg_sm, tadpole{1,t}.multimax_peakavg_sm)
    xlabel('unisensory avg peak response')
    ylabel('multisensory avg peak response')
    title(sprintf('tadpole %d peak response', t))
    fig_filename=sprintf('peak_univmulti_tadpole%d_smoothed', t);
    saveas(gcf,fig_filename,'png');
end

% scatter all together
univals = [];
multivals = [];
for t = 1:length(tadpole)
    univals = [univals tadpole{1,t}.unimax_peakavg_sm]
    multivals = [multivals tadpole{1,t}.multimax_peakavg_sm]
end
scatter(univals, multivals)
xlabel('unisensory avg peak response')
ylabel('multisensory avg peak response')
title('All tadpoles peak response')
fig_filename='peak_univmulti_all_smoothed';
saveas(gcf,fig_filename,'png')

%% Peak location histogram

% multi
multi_peaklocall = [];
for t = 1:length(tadpole)
    multi_peaklocall = [multi_peaklocall tadpole{1,t}.peakloc_avg_sm(1,:)];
end
figure;
hist(multi_peaklocall,40)

uni_peaklocall = [];
for t = 1:length(tadpole)
    for i = 1:size(tadpole{1,t}.peakloc_avg_sm,2)
        if tadpole{1,t}.unimax_stimtype_sm(1,i) == 1 %visual
            uni_peaklocall = [uni_peaklocall tadpole{1,t}.peakloc_avg_sm(2,i)];
        elseif tadpole{1,t}.unimax_stimtype_sm(1,i) == 2
            uni_peaklocall = [uni_peaklocall tadpole{1,t}.peakloc_avg_sm(3,i)];
        end
    end
end
figure;
hist(uni_peaklocall,40)

difference = multi_peaklocall - uni_peaklocall;
figure;
hist(difference,50)
[h, p] = kstest(difference)
% h=1, so difference is not normally distributed

%% For every cell, percent response for each category

%% For each cell, peak response
for t = 1:length(tadpole)
    for i = 1:size(tadpole{1,t}.peak_avg,2)
        tmpmax = max(tadpole{1,t}.peak_avg(1,i)-tadpole{1,t}.peak_avg(4,i));
        spontdiff(i) = tmpmax-tadpole{1,t}.peak_avg(4,i);
    end
    tadpole{1,t}.spontdiff = spontdiff;
    clear('spontdiff')
end

% how many cells have positive spontdiff (indicating that at least one mean stimulus is
% greater than mean spontaneous activity)

for t = 1:length(tadpole)
    pos_spontdiff(t) = length(find(tadpole{1,t}.spontdiff > 0));
    neg_spontdiff(t) = length(find(tadpole{1,t}.spontdiff < 0));
end
total_pos = sum(pos_spontdiff)
total_neg = sum(neg_spontdiff)
prop_posspontdiff = pos_spontdiff ./ (neg_spontdiff + pos_spontdiff)

%% Proximity 
% scatter plot and linear fit of euclid distance vs peak avg loc

