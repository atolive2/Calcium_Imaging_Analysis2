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
            if multi > 0
                tadcluster{1,t}.MSEnh_numresponses(r) = 1;
            else
                tadcluster{1,t}.MSEnh_numresponses(r) = 0; 
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

%% Get max Rvals and lag times so can compare

% Make df_f0 separated by trial type (only use 1-4 here to eliminate
% variation from stim strength)
tmp = [];
for t = 1:length(tadcluster)
    for s = 1:4
        trials_touse = find(tadcluster{1,t}.stimorder == s)
       %tmp = [];
        
            for j = 1:size(tadcluster{1,t}.df_f0,1) % over each ROI
                for i = 1:length(trials_touse) %over all trials of 1 stim type
                    tmp = [tmp; tadcluster{1,t}.df_f0{j, trials_touse(i)}];
                end
                tadcluster{1,t}.dff0_bystimtype{j, s} = tmp;
                tmp = [];
            end
    end
end

% first get a matrix ROI x dff0 for all trials of a given stimtype
for t = 1:length(tadcluster)
    for i = 1:size(tadcluster{1,t}.dff0_bystimtype,1)
        if sum(size(tadcluster{1,t}.dff0_bystimtype{i,1})) > 0 %elimnate tads with no high/high multisensory trials
        tadcluster{1,t}.dff0_multi(i,:) = tadcluster{1,t}.dff0_bystimtype{i,1}';
        tadcluster{1,t}.dff0_vis(i,:) = tadcluster{1,t}.dff0_bystimtype{i,2}';
        tadcluster{1,t}.dff0_mech(i,:) = tadcluster{1,t}.dff0_bystimtype{i,3}';
        tadcluster{1,t}.dff0_none(i,:) = tadcluster{1,t}.dff0_bystimtype{i,4}';
        end
    end
end

% Calculate correlation coefficients by stim type

for t = 1:length(tadcluster)
    if sum(size(tadcluster{1,t}.dff0_bystimtype{1,1})) > 0
        set_lag = length(tadcluster{1,t}.df_f0{1,1});
        [tadcluster{1,t}.respROIdff0_R_MS, tadcluster{1,t}.respROIdff0_lag_MS] = xcorr(tadcluster{1,t}.dff0_multi(tadcluster{1,t}.resp_ROIs,:)', set_lag, 'coeff');
        [tadcluster{1,t}.respROIdff0_R_V, tadcluster{1,t}.respROIdff0_lag_V] = xcorr(tadcluster{1,t}.dff0_vis(tadcluster{1,t}.resp_ROIs,:)', set_lag, 'coeff');
        [tadcluster{1,t}.respROIdff0_R_M, tadcluster{1,t}.respROIdff0_lag_M] = xcorr(tadcluster{1,t}.dff0_mech(tadcluster{1,t}.resp_ROIs,:)', set_lag, 'coeff');
        [tadcluster{1,t}.respROIdff0_R_N, tadcluster{1,t}.respROIdff0_lag_N] = xcorr(tadcluster{1,t}.dff0_none(tadcluster{1,t}.resp_ROIs,:)', set_lag, 'coeff');
    end
end

% what is the max correlation between 2 ROIs?
% this is the max of each column, and all columns that are an
% autocorrelation will have max = 1
for t = 1:length(tadcluster)
    if sum(size(tadcluster{1,t}.dff0_bystimtype{1,1})) > 0
    [tadcluster{1,t}.respROIdff0_maxR_MS(1,:), tadcluster{1,t}.respROIdff0_maxR_MS(2,:)] = max(tadcluster{1,t}.respROIdff0_R_MS);
    [tadcluster{1,t}.respROIdff0_maxR_V(1,:), tadcluster{1,t}.respROIdff0_maxR_V(2,:)] = max(tadcluster{1,t}.respROIdff0_R_V);
    [tadcluster{1,t}.respROIdff0_maxR_M(1,:), tadcluster{1,t}.respROIdff0_maxR_M(2,:)] = max(tadcluster{1,t}.respROIdff0_R_M);
    [tadcluster{1,t}.respROIdff0_maxR_N(1,:), tadcluster{1,t}.respROIdff0_maxR_N(2,:)] = max(tadcluster{1,t}.respROIdff0_R_N);
    end
end

%reshape the array to be size(length(resp_ROIs)) e.g. make it a square
%again for plotting purposes. 
%for reference, xcorr arranges the cols as 1-1, 1-2, 1-3 ... 1-n, 2-1, 2-2,
%etc
for t = 1:length(tadcluster)
    if sum(size(tadcluster{1,t}.dff0_bystimtype{1,1})) > 0
        len = length(tadcluster{1,t}.resp_ROIs);
        for r = 1:len
            for c = 1:len
                idx = (r-1)*len + c;
                % max R value
                tadcluster{1,t}.respROIdff0_maxR_sq_MS(r,c) = tadcluster{1,t}.respROIdff0_maxR_MS(1,idx);
                tadcluster{1,t}.respROIdff0_maxR_sq_V(r,c) = tadcluster{1,t}.respROIdff0_maxR_V(1,idx);
                tadcluster{1,t}.respROIdff0_maxR_sq_M(r,c) = tadcluster{1,t}.respROIdff0_maxR_M(1,idx);
                tadcluster{1,t}.respROIdff0_maxR_sq_N(r,c) = tadcluster{1,t}.respROIdff0_maxR_N(1,idx);
                % lag time of max R val
                tadcluster{1,t}.respROIdff0_maxRlag_sq_MS(r,c) = tadcluster{1,t}.respROIdff0_lag_MS(1,tadcluster{1,t}.respROIdff0_maxR_MS(2,idx));
                tadcluster{1,t}.respROIdff0_maxRlag_sq_V(r,c) = tadcluster{1,t}.respROIdff0_lag_V(1,tadcluster{1,t}.respROIdff0_maxR_V(2,idx));
                tadcluster{1,t}.respROIdff0_maxRlag_sq_M(r,c) = tadcluster{1,t}.respROIdff0_lag_M(1,tadcluster{1,t}.respROIdff0_maxR_M(2,idx));
                tadcluster{1,t}.respROIdff0_maxRlag_sq_N(r,c) = tadcluster{1,t}.respROIdff0_lag_N(1,tadcluster{1,t}.respROIdff0_maxR_N(2,idx));
            end
        end
    end
end

% Go run lines 465-625 of correlations_newdata to find the highly
% correlated ROIs

%% ID highly correlated ROIs
%%%%%%%%%%%%%%%% PROBLEM. HIGH CORR CELLS NOT APPEARING IN TADS 1-15. 

for t = 1:length(tadcluster)
    if isfield(tadcluster{1,t}, 'correlated_ROIs_dff0_MS_common_AROI')
        high_corrs = [];
        sprintf('tad %d has high corrs', t)
        for i = 1:length(tadcluster{1,t}.correlated_ROIs_dff0_MS_common_AROI)
            high_corrs = [high_corrs; tadcluster{1,t}.correlated_ROIs_dff0_MS_common_AROI{i}];
        end
        for i = 1:length(tadcluster{1,t}.resp_ROIs)   
            if length(find(high_corrs == tadcluster{1,t}.resp_ROIs(i))) > 0
                tadcluster{1,t}.ishighcorr_MS(i) = 1;
            else 
                tadcluster{1,t}.ishighcorr_MS(i) = 0;
            end
        end
    else
        tadcluster{1,t}.ishighcorr_MS = zeros(length(tadcluster{1,t}.resp_ROIs));
    end
end

for t = 1:length(tadcluster)
    if isfield(tadcluster{1,t}, 'ishighcorr_MS')
        sprintf('%d', t)
    end
end

%% Put data into 1 matrix for clustering
% only use responding ROIs

% to add together multiple tadclusters - rename vars first
%tadcluster = [tadcluster2, tadcluster1];


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
        cluster_data_respROI(counter, 22) = tadcluster{1,t}.ishighcorr_MS(rr); %IS the ROI "highly correlated" from xcorr multi
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
hist(SigP(:,2), 20)
xlim([1, 20])
xlabel('var num')
ylabel('num correlated')
title('Correlated Variables across all cells')
fig_filename = 'hist of correlated input variables new'
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

%% Plot MDS coordinates colored by primary modality
% use enoughROItads to calculate primary modality as the maximum average
% peak
for i = 1:length(enoughROItads)
    [PrimaryMod{i}(:,1), PrimaryMod{i}(:,2)] = max(enoughROItads{i}(:,3:6),[],2);
end

% generate color list based on primary modality
for i = 1:length(PrimaryMod)
    for j = 1:size(PrimaryMod{i})
    if PrimaryMod{i}(j,2) == 1 %multi
        Colors{i}(j,:) = [1 0 1]; %purple
    elseif PrimaryMod{i}(j,2) == 2 %vis
        Colors{i}(j,:) = [0 0 1]; %blue
    elseif PrimaryMod{i}(j,2) == 3 %mech
        Colors{i}(j,:) = [1 0 0]; %red
    elseif PrimaryMod{i}(j,2) == 4 %no stim
        Colors{i}(j,:) = [0.5 0.5 0.5]; %gray
    end
    end
end

%plot Y with colors from Colors
for i = 1:length(Y)
    figure;
    scatter(Y_B{i}(:,1), Y_B{i}(:,2), 50, Colors{i}, 'filled')
    title(sprintf('tad %d (%d) 2D MDS respROIs colored by Primary Modality', i, enoughROItads{i}(1,1)))
    xlabel('MDS dim 1')
    ylabel('MDS dim 2')
    annotation('textbox', 'Position', [0.7 0.75 .1 .1], 'String', ['Multi'], 'Color', 'm', 'LineStyle', 'none' );
    annotation('textbox', 'Position', [0.7 0.7 .1 .1], 'String', ['Vis'], 'Color', 'b', 'LineStyle', 'none' );
    annotation('textbox', 'Position', [0.7 0.65 .1 .1], 'String', ['Mech'], 'Color', 'r', 'LineStyle', 'none' );
    fig_filename = sprintf('tad %d (%d) 2D MDS respROIs colored by Primary Modality', i, enoughROItads{i}(1,1))
    saveas(gcf, fig_filename, 'png')
    close;
end

%% Plot relationship between topographical distance and MDS distance 
% all cells on 1 plot, color by experiment

% get MDS distance between ROIs
for i = 1:length(Y)
     Y_Distance{i} = pdist(Y{i})
end

%get topographical distance between ROIs
for i = 1:length(Y)
    respROI_dist{i} = tadcluster{1, enoughROItads{i}(1,1)}.ROIcenters(enoughROItads{i}(:,2), :)
    T_Distance{i} = pdist(respROI_dist{i});
end

%normalize data by Z score
for i = 1:length(Y_Distance)
    Y_distance_Z{i} = zscore(Y_Distance{i});
    T_distance_Z{i} = zscore(T_Distance{i});
end

% make trendline
% assemble data into 1 long vector
Y_distanceZ_all = [];
T_distanceZ_all = [];
for i = 1:length(Y_distance_Z)
    Y_distanceZ_all = [Y_distanceZ_all Y_distance_Z{i}];
    T_distanceZ_all = [T_distanceZ_all T_distance_Z{i}];
end

%generate fit data points
my_poly=polyfit(Y_distanceZ_all, T_distanceZ_all, 1);
X2= min(Y_distanceZ_all):0.1:max(Y_distanceZ_all); % X data range 
Y2=polyval(my_poly,X2);

plot(X2,Y2);


% build scatterplot
cmap = colormap(jet(length(T_Distance)));
figure;
hold on
for i = 1:length(Y_Distance)
    scatter(Y_distance_Z{i}, T_distance_Z{i}, 50, cmap(i,:), 'filled')
end
plot(X2,Y2, 'k', 'LineWidth', 5);
hold off
title('MDS Distance vs Topographical Distance by Exp')
xlabel('MDS distance (Z-score)')
ylabel('Topographical Distance (Z-score)')
annotation('textbox', 'Position', [0.7 0.75 .1 .1], 'String', ['Y = 0.3X + ~0'], 'Color', 'k', 'LineStyle', 'none' );
fig_filename = 'MDS distance vs topographical distance all tads'
saveas(gcf, fig_filename, 'png')

%Plot each experiment separately
for i = 1:length(Y_distance_Z)
    %generate fit data points
    my_poly=polyfit(Y_distance_Z{i}, T_distance_Z{i}, 1);
    X2= min(Y_distance_Z{i}):0.1:max(Y_distance_Z{i}); % X data range 
    Y2=polyval(my_poly,X2);
    figure;
    hold on
    scatter(Y_distance_Z{i}, T_distance_Z{i}, 50, 'filled')
    plot(X2,Y2,'k', 'LineWidth', 5);
    annotation('textbox', 'Position', [0.5 0.8 .1 .1], 'String', ['Y = ', num2str(my_poly(1,1)),'X + ', num2str(my_poly(1,2))], 'Color', 'k', 'LineStyle', 'none' );
    title(sprintf('tad %d (%d) MDS Distance vs Topographical Distance', i, enoughROItads{i}(1,1)))
    xlabel('MDS distance (Z-score)')
    ylabel('Topographical Distance (Z-score)')
    fig_filename = sprintf('tad %d (%d) MDS distance vs topographical distance', i, enoughROItads{i}(1,1))
    saveas(gcf, fig_filename, 'png')
    close;
    clear('mypoly', 'X2', 'Y2')
end


%% Plot MDS colored by Multisensory enhancement
%MSEnh by peak is stored in enoughROItads{i}(:,17)

% generate color list based on MSEnh
for i = 1:length(enoughROItads)
    MSEnh = enoughROItads{i}(:,17)
    [MSEnh_sort{i}(:,1), MSEnh_sort{i}(:,2)] = sortrows(MSEnh) %sort from smallest to largest MSEnh 
    cmap = colormap(jet(size(enoughROItads{i},1))); %get a list of colors
    MSEnh_sort{i}(:,3:5) = cmap %add the cmap to the MSEnh vals in order from smallest to largest
    MSEnh_resort{i} = sortrows(MSEnh_sort{i}, 2) %sort by the original sortrows index to get back to MDS data order
end

% Make scatter plot to plot MDS coordinates colored by MSEnh
for i = 1:length(Y)
    figure;
    scatter(Y{i}(:,1), Y_B{i}(:,2), 50, MSEnh_resort{i}(:,3:5), 'filled')
    title(sprintf('tad %d (%d) 2D MDS respROIs colored by MSEnh', i, enoughROItads{i}(1,1)))
    xlabel('MDS dim 1')
    ylabel('MDS dim 2')
    fig_filename = sprintf('tad %d (%d) 2D MDS respROIs colored by MSEnh', i, enoughROItads{i}(1,1))
    saveas(gcf, fig_filename, 'png')
    close;
end
    
%% Plot MDS colored by mean onset time
%Mean onset time is stored in enoughROItads{i}(:,11:14)
%edited to be  by multi, then vis, then mech, individually

% generate color list based on Mean onset time
for i = 1:length(enoughROItads)
    MSEnh = enoughROItads{i}(:,13)
    [MSEnh_sort{i}(:,1), MSEnh_sort{i}(:,2)] = sortrows(MSEnh) %sort from smallest to largest MSEnh 
    cmap = colormap(jet(size(enoughROItads{i},1))); %get a list of colors
    MSEnh_sort{i}(:,3:5) = cmap %add the cmap to the MSEnh vals in order from smallest to largest
    MSEnh_resort{i} = sortrows(MSEnh_sort{i}, 2) %sort by the original sortrows index to get back to MDS data order
end

% Make scatter plot to plot MDS coordinates colored by MSEnh
for i = 1:length(Y)
    figure;
    scatter(Y{i}(:,1), Y_B{i}(:,2), 50, MSEnh_resort{i}(:,3:5), 'filled')
    title(sprintf('tad %d (%d) 2D MDS respROIs colored by Mean Onset Time (Mech)', i, enoughROItads{i}(1,1)))
    xlabel('MDS dim 1')
    ylabel('MDS dim 2')
    fig_filename = sprintf('tad %d (%d) 2D MDS respROIs colored by Mean Onset Time (Mech)', i, enoughROItads{i}(1,1))
    saveas(gcf, fig_filename, 'png')
    %close
end
