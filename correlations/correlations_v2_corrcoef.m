%% Correlations Version 2: calculate correlation coefficients based on random data about ROIs using corrcoef
% R = corrcoef(A) returns the matrix of correlation coefficients for A, 
% where the columns of A represent random variables and the rows represent observations.
% https://www.mathworks.com/help/matlab/ref/corrcoef.html

% all of this is copied from correlations_bymodality and
% cluster_analysis_v1 but using random variable data instead of time series

% this code starts with tadcluster_data_20170814, with all fields after
% resp_ROIs deleted.

%% What fields do I want to include? 
% basically, let's put everything in and see what happens. 

%% Copied from cluster_analysis_v2
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

%% Use cluster_data_respROI to calculate correlation coefficients of the input variables 

% this produces symmetric matrix of correlation coeeficients nad the matrix
% of Pvals for testing the hypothesis that thereis no relationship between
% the observed phenomenon.
[R, P] = corrcoef(cluster_data_respROI(:,3:end), 'rows', 'pairwise')

% are there any significant relationships between variables?
[SigP(:,1), SigP(:,2)] = find(P < 0.05) 

% remove "duplicates" from SigP (that is, eliminate 19-1 because 1-19
% already exists)
SigP_sm = [];
for i = 1:size(SigP,1)
    if SigP(i,1) > SigP(i,2)
        SigP_sm = [SigP_sm; SigP(i,:)] 
    end
end

%[a,b] = hist(SigP_sm(:,2), unique(SigP_sm(:,2)))
hist(SigP(:,2), 19)
xlim([1, 19])
fig_filename = 'hist of correlated input variables'
saveas(gcf, fig_filename, 'png')
%so the variables not correlated with any other variables are
%8 (peak_stddev_bymod Vis), 18 (MSEnh_numresponses) and 19
%(unibias_num_responses) 
%hunch is that MSEnh and unibias have too many NaNs to be useful
% some vars are correlated with multiple things
% so now what? I think I try to eliminate some variables and use only
% variables that differentiate the cells in different ways a cluster analysis

%% Let's try: Classical Multidimensional Scaling
% this is a method of nonlinear dimensionality reduction that visualizes
% the level of similarity of individual cases in a dataset. 
% https://www.mathworks.com/help/stats/multidimensional-scaling.html
% https://www.mathworks.com/help/stats/cmdscale.html
% this paper does it in monkeys: https://www.ncbi.nlm.nih.gov/pubmed/25728571
% basics: https://en.wikipedia.org/wiki/Multidimensional_scaling

% first, create the dissimilarity matrix
% You can specify D as either a full dissimilarity matrix, or in upper 
% triangle vector form such as is output by pdist. 
% use pdist on portions of pdist that correspond to tadpoles with enough
% ROIs
% https://www.mathworks.com/help/stats/pdist.html

% determine the number of ROIs per tad in cluster_data_respROI
tads = unique(cluster_data_respROI(:,1))
for t = 1:length(tads)
    roi_idx{t} = find(cluster_data_respROI(:,1) == tads(t))
end

%eliminate tads with less than 18 ROIs
counter=1
for i = 1:length(roi_idx)
    if length(roi_idx{i}) >= 18
        enoughROItads{counter} = cluster_data_respROI(roi_idx{i},:);
        counter = counter + 1;
    end
end

% create distance matrix for each tad with 18+ ROIs
% nan is a problem here. Don't see obvious solution so removing last 2
% columns (MSEnh_numrepsonses and unibias_numresponses) for now. 
% default for pdist is euclidean distance.
for i = 1:length(enoughROItads)
    D{i} = pdist(enoughROItads{i}(:,3:19))
end

% feed distance matrix into cmdscale
for i = 1:length(D)
    Y{i} = cmdscale(D{i},2)
end

% sort the rows of Y to allow coloring by location to be useful
for i = 1:length(Y)
    [Y_B{i}, roi_order{i}] = sortrows(Y{i});
    
end

% plot the MDS for each tadpole
for i = 1:length(Y_B)
    figure;
    c = linspace(1,10, size(Y_B{i},1));
    scatter(Y_B{i}(:,1), Y_B{i}(:,2), 50, c, 'filled')
    clear('c')
    title(sprintf('tad %d (%d) 2D MDS respROIs', i, enoughROItads{i}(1,1)))
    xlabel('MDS dim 1')
    ylabel('MDS dim 2')
    %fig_filename = sprintf('tad %d (%d) 2D MDS respROIs', i, enoughROItads{i}(1,1))
    %saveas(gcf, fig_filename, 'png')
    %close;
end

% plot the ROIs tectum shaped and colored same as MDS
% align real ROI locations with order of plotting from MDS
for i = 1:length(roi_order)
    roi_labels = enoughROItads{i}(:,2);
    real_rois = roi_labels(roi_order{i});
    roi_order{i}(:,2:3) = tadcluster{1, enoughROItads{i}(1,1)}.ROIcenters(real_rois, :);
    clear('roi_labels', 'real_rois')
end

% make scatter plot
for i = 1:length(roi_order)
    t = enoughROItads{i}(1,1);
    figure;
    hold on
    scatter(tadcluster{1,t}.ROIcenters(:,1), tadcluster{1,t}.ROIcenters(:,2), '+', 'MarkerEdgeColor', 'k') %plot all ROIs
    c = linspace(1,10, size(Y_B{i},1)); % create same color vector as MDS
    scatter(roi_order{i}(:,2), roi_order{i}(:,3), 50, c, 'filled') %plot respROIs in colors of MDS
    hold off
    title(sprintf('tad %d (%d) tectum by respROIs MDS color', i, enoughROItads{i}(1,1)))
end

