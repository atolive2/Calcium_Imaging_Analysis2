%% Cluster Analysis version 2 [by modality]
% Use specific input variables to make clusters 
% Normalize to quartiles/median to decrease effects of outliers. 

% begin with tadcluster_data_20170814.mat (this is just tadcluster), and
% deleted all fields after resp_ROIs.

% calculate variables by cell. Include responding ROIs only. 
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

%%%%%%% for me to clear vars when I screw up
for t = 1:length(tadcluster)
tadcluster{1,t}.peak_stddev_bymod = [];
tadcluster{1,t}.peak_stddev_MS = [];
end
%%%%%%%%%

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
            tadcluster{1,t}.MSEnh_numresponses(r) = NaN;
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
            tadcluster{1,t}.unibias_numresponses(r) = NaN;
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

%% Run cluster analysis

T = clusterdata(cluster_data_respROI(:, 3:19),20)
% cutoff: 1 = 217 clusters, 1.1 = 213, 1.2:1.9 = 1 cluster, 2 = 2 clusters,
% 20 = 20 clusters.
% What does this mean???????????

% PCA modelled from https://www.mathworks.com/products/demos/machine-learning/cluster_genes/cluster_genes.html
[~,score,~,~,explainedVar] = pca(cluster_data_respROI(:, 3:19));
bar(explainedVar)
title('Explained Variance: More than 90% explained by first three principal components')
ylabel('PC')

% Retain first 3 principal components
tadPC = score(:,1:3);

figure;
[clusters, centroid] = kmeans(tadPC,6);
scatter3(tadPC(:,1),tadPC(:,2),tadPC(:,3),clusters)
legend('location','southeast')
xlabel('First Principal Component');
ylabel('Second Principal Component');
zlabel('Third Principal Component');
title('Principal Component Scatter Plot with Colored Clusters');

% neural network clustering
%input data can't have NaN
cluster_data_NNcluster = cluster_data_respROI(:, 3:19);
%% do clusters map on to:
% - primary modality (across tadpoles)
% - physical space/proximity to similar cells (by tadpole)

% what input properties drive cluster formation?
% maybe a PCA would be useful?