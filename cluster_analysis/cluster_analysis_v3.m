%% Cluster Analysis Version 3: By modality, with range normalization by IQR
% this will hopefully remove the impact of outliers on normalization (a
% problem with regular range normalization and Z scoring)

%started with tadcluster_data_20170814_clean.mat which is all the basic ROI
%info

%% First collect all relevant parameters into 1 matrix (same as cluster_analysis_v2)

% Carlos suggests: avg peak, jitter, mean squared error, hindbrain bias, peaks by modality
% overall strategy: calculate for all cells, then extract just the cells I
% want to cluster at the end.

%% peaks by modality - mean and stddev of each stim type
% avg_peak is in tadcluster{1,t}.peak_avg

for t = 1:length(tadcluster)
    tadcluster{1,t}.stimmask = get_stimmask(tadcluster{1,t}.stimorder)
end
for t = 1:length(tadcluster)
    tadcluster{1,t}.peak_stddev_bymod = std_by_stimtype(tadcluster{1,t}.peak_bytrial, tadcluster{1,t}.stimmask)
end

%% jitter = variation in onset time (separate by modality)
% copy code from analyze_onset_time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%add stim onset time to tadcluster%%%%
% exps 1-9 have onset at 0.5s 
% exps 10+ have onset at 2s
for t = 1:length(tadcluster)
    if tadcluster{1,t}.expnum <=9
        tadcluster{1,t}.stim_onset = 0.5;
    elseif tadcluster{1,t}.expnum >9
        tadcluster{1,t}.stim_onset = 2;
    else
        fprintf('error exp %d, t')
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get time to 50% of the peak (= response onset time)
% same as used in analyze_onset_time
for t = 1:length(tadcluster)
    for r = 1:size(tadcluster{1,t}.df_f0,1) %over each ROI
        for tt = 1:size(tadcluster{1,t}.df_f0,2) %over each trial
            half_peak = tadcluster{1,t}.peak_bytrial(r, tt) / 2;
            tmp_onset = find((tadcluster{1,t}.df_f0{r, tt} > half_peak) , 3);
            % get trial length (this takes care of exps with not 160 frames)
            trial_length = length(tadcluster{1,t}.df_f0{r, tt});
            %find the proper onset time (e.g. eliminate errors of 1 or 2)
            if length(tmp_onset) < 3 %this sets onset time to 0 for trials with no response
                onset_time(r, tt) = 0;
            else
                if tmp_onset(1) == (1 || 2)
                    if tmp_onset(2) == 2 
                        onset_time(r, tt) = tmp_onset(3) / (trial_length/7);
                    else
                        onset_time(r, tt) = tmp_onset(2) / (trial_length/7);
                    end
                else
                    onset_time(r, tt) = tmp_onset(1) / (trial_length/7);
                end
            end
        end
    end
    % put onset_time into tadcluster
    tadcluster{1,t}.onset_time = onset_time;    
end

% Get mean and sd of onset_time by ROI (this includes 0s - no response trials)
for t = 1:length(tadcluster)
    tadcluster{1,t}.onset_time_mean = mean_by_stimtype2(tadcluster{1,t}.onset_time, tadcluster{1,t}.stimmask);
    tadcluster{1,t}.onset_time_stdev = std_by_stimtype(tadcluster{1,t}.onset_time, tadcluster{1,t}.stimmask);
end

%% multisensory enhancement (recalculate by num responses)

% use tadcluster{1,t}.MSEnh_peak 

% number of responses
for  t = 1:length(tadcluster)
    for r = 1:size(tadcluster{1,t}.boolean_response, 1)
        multi = sum(tadcluster{1,t}.boolean_response(r, tadcluster{1,t}.stimmask(:,1))) / sum(tadcluster{1,t}.stimmask(:,1));
        vis = sum(tadcluster{1,t}.boolean_response(r, tadcluster{1,t}.stimmask(:,2))) / sum(tadcluster{1,t}.stimmask(:,2));
        mech = sum(tadcluster{1,t}.boolean_response(r, tadcluster{1,t}.stimmask(:,3))) / sum(tadcluster{1,t}.stimmask(:,3));
        if (vis + mech) == 0 
            if multi > 0
                tadcluster{1,t}.MSEnh_numresponses(r) = 1;
            else
                tadcluster{1,t}.MSEnh_numresponses(r) = NaN; 
            end
        else
            tadcluster{1,t}.MSEnh_numresponses(r) = multi / (vis + mech);
        end
    end
end

%% unisensory bias = ratio of vis responses to mech responses

% get number of responses to vis and mech, divide vis by mech
for  t = 1:length(tadcluster)
    for r = 1:size(tadcluster{1,t}.boolean_response, 1)
        vis = sum(tadcluster{1,t}.boolean_response(r, tadcluster{1,t}.stimmask(:,2))) / sum(tadcluster{1,t}.stimmask(:,2));
        mech = sum(tadcluster{1,t}.boolean_response(r, tadcluster{1,t}.stimmask(:,3))) / sum(tadcluster{1,t}.stimmask(:,3));
        if (vis + mech) == 0 %doesn't respond to either
            tadcluster{1,t}.unibias_numresponses(r) = 0;
        else
            tadcluster{1,t}.unibias_numresponses(r) = vis/(vis+mech);
        end
    end
end

%% Put data into 1 matrix for clustering
% only use responding ROIs
counter = 1
for t = 1:length(tadcluster)
    for rr = 1:length(tadcluster{1,t}.resp_ROIs) 
        roi = tadcluster{1,t}.resp_ROIs(rr);
        cluster_data_respROI(counter,1) = t; %tad num
        cluster_data_respROI(counter,2) = roi; %roi num
        cluster_data_respROI(counter,3:6) = tadcluster{1,t}.peak_avg(1:4,roi)'; %avg peak over high/high stimuli
        cluster_data_respROI(counter,7:10) = tadcluster{1,t}.peak_stddev_bymod(1:4,roi)'; %stdev of peak over high/high stimuli
        cluster_data_respROI(counter,11:14) = tadcluster{1,t}.onset_time_mean(1:4,roi)'; %avg onset time over high/high stimuli [Jitter]
        cluster_data_respROI(counter,15:18) = tadcluster{1,t}.onset_time_stdev(1:4,roi)'; %stdev onset time over high/high stimuli [Jitter]
        cluster_data_respROI(counter,19) = tadcluster{1,t}.MSenh_peak(roi); %Multisensory enhancement based on average peak over high/high stimuli
        cluster_data_respROI(counter,20) = tadcluster{1,t}.MSEnh_numresponses(roi); %Multisensory enhancement based on number of responses to high/high stimuli
        cluster_data_respROI(counter,21) = tadcluster{1,t}.unibias_numresponses(roi); %unisensory bias based on number of responses to high/high stimuli
        counter = counter + 1 ;
    end
end

unique(cluster_data_respROI(:,1)) %all exps represented except for 13

%% Normalize each variable by IQR
% set median = 0 and 25% as -1 and 75% as 1

for i = 1:size(cluster_data_respROI,2) %for each variable
    rng = iqr(cluster_data_respROI(:,i))
    cluster_data_norm(:,i) = 
end

cluster_data_iqr = iqr(cluster_data_respROI,1)
cluster_data_med = median(cluster_data_respROI,1)
for i = 1:length(cluster_data_iqr)
    cluster_data_norm(:,i) = cluster_data_respROI(:,i) / cluster_data_iqr(i);
    iqr_check(i) = iqr(cluster_data_norm(:,i));
end
cluster_data_med = median(cluster_data_norm,1);

%plot raw data vs IQR normed data (check and ensure no data integrity lost)
for i = 1:size(cluster_data_respROI,2)
    figure;
    scatter(cluster_data_respROI(:,i), cluster_data_norm(:,i),50, 'filled')
end
%all look good. Outliers still look like outliers, everything else is
%clustered in a smaller range. 

%% Do cluster analysis using normalized data

% version 1: heirarchacal clustering with dendrograms
% this is usually used on genes: https://www.mathworks.com/help/bioinfo/ref/clustergram.html

%NaN is not allowed by this function
CGobj = clustergram(cluster_data_norm(:,[1:19,21]))
% result is kind of useless because it's for every responding ROi across
% all exps

%try clustergram separately on each tad
% determine the number of ROIs per tad in cluster_data_respROI
tads = unique(cluster_data_norm(:,1))
for t = 1:length(tads)
    roi_idx{t} = find(cluster_data_norm(:,1) == tads(t))
end
%eliminate tads with less than 18 ROIs
counter=1
for i = 1:length(roi_idx)
    if length(roi_idx{i}) >= 18
        enoughROItads{counter} = cluster_data_norm(roi_idx{i},:);
        counter = counter + 1;
    end
end

labels = {'peak MS', 'peak V', 'peak M', 'peak NS', 'peak SD MS', 'peak SD V', 'peak SD M', 'peak SD NS', 'onset time MS', 'onset time V', 'onset time M', 'onset time NS', ...
    'onset time SD MS', 'onset time SD V', 'onset time SD M', 'onset time SD NS', ...
    'MSEnh peak', 'Uni bias num resp'}
for i = 1:length(enoughROItads)
    CGobj = clustergram(enoughROItads{i}(:,[3:19,21]), 'columnlabels', labels)
    Cluster_obj{i} = CGobj
    clear('CGobj')
end

%%%%% This didn't really work. 

%% Let's plot all extracted vals by ROI for all responding ROIs

labels = {'peak MS', 'peak V', 'peak M', 'peak NS', 'peak SD MS', 'peak SD V', 'peak SD M', 'peak SD NS', 'onset time MS', 'onset time V', 'onset time M', 'onset time NS', ...
    'onset time SD MS', 'onset time SD V', 'onset time SD M', 'onset time SD NS', ...
    'MSEnh peak', 'MSEnh num resp', 'Uni bias num resp'}

for i = 3:size(cluster_data_respROI, 2) %over each extracted parameter
    for j = 3:size(cluster_data_respROI, 2) %over each extracted parameter again
        if i ~= j
            rois_touse = [];
            for k = 1:size(cluster_data_respROI, 1)
                if (max(cluster_data_respROI(k,[i,j])) < 3) && (min(cluster_data_respROI(k,[i,j])) > -3)
                    rois_touse = [rois_touse, k]; %3 zscores above or below 0 = outlier 
                end
            end
        p = polyfit(cluster_data_respROI(rois_touse,i), cluster_data_respROI(rois_touse,j),1)
        X2 = min(cluster_data_respROI(rois_touse,i)):0.1:max(cluster_data_respROI(rois_touse,i)); % X data range 
        Y2 = polyval(p,X2);
        figure;
        hold on
        plot(cluster_data_respROI(rois_touse,i), cluster_data_respROI(rois_touse,j), 'ok')
        plot(X2,Y2, 'b', 'LineWidth', 1);
        xlabel(labels{i-2})
        ylabel(labels{j-2})
        title(sprintf('No outliers All Resp ROI %s vs %s', labels{i-2}, labels{j-2}))
        fig_filename = sprintf('No outliers All Resp ROI %s vs %s', labels{i-2}, labels{j-2})
        saveas(gcf, fig_filename, 'png')
        close;
        clear('p', 'X2', 'Y2', 'fig_filename');
        end
    end
end

% are any significant?
for i = 3:size(cluster_data_respROI, 2) %over each extracted parameter
    for j = 3:size(cluster_data_respROI, 2)
        rois_touse = [];
        for k = 1:size(cluster_data_respROI, 1)
            if (max(cluster_data_respROI(k,[i,j])) < 3) && (min(cluster_data_respROI(k,[i,j])) > -3)
                rois_touse = [rois_touse, k] %3 zscores above or below 0 = outlier 
            end
        end
        [R_respROI(i-2,j-2), P_respROI(i-2,j-2)] = corrcoef(cluster_data_respROI(rois_touse,i), cluster_data_respROI(rois_touse,j))
    end
end

[R_respROI_Woutliers, P_respROI_Woutliers] = corrcoef(cluster_data_respROI, 'rows','pairwise')
SigP = logical(P_respROI_Woutliers < 0.05)
highR = logical(R_respROI_Woutliers > 0.5)
intersection = SigP & highR
[Sig_combo(:,1), Sig_combo(:,2)] = find(intersection)
size(Sig_combo) %yields 12 -> 6 pairs of vars with significant P vals and high R vals
for i = 1:size(Sig_combo, 1)
    for j = 1:size(Sig_combo, 2)
        Sig_combo_names{i,j} = labels{Sig_combo(i,j)-2}
    end
end
for i = 1:size(Sig_combo, 1)
        Sig_combo_vals(i,1) = R_respROI_Woutliers(Sig_combo(i,1), Sig_combo(i,2))
        Sig_combo_vals(i,2) =  P_respROI_Woutliers(Sig_combo(i,1), Sig_combo(i,2))
end

%% A different cluster analysis using normalized data
