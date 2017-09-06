%% Examples for retreat talk 
% tadpole 12 - 110 ROIs, clear stripy-ness. 

% open file 
load('tadcluster_analysis_xcorr_20170818.mat')

%% first, replot correlation maxR vals with the same color axis values



%% Second, find cells that have a high correlation with at least 25% of other cells
t = 12
% find ROIs with at least 0.25*number of ROIs highly correlated
high_MS = find(tadcluster{1,t}.highcorr_numROIs_MS > (0.25*length(tadcluster{1,t}.highcorr_numROIs_MS)))
high_V = find(tadcluster{1,t}.highcorr_numROIs_V > (0.25*length(tadcluster{1,t}.highcorr_numROIs_V)))
high_M = find(tadcluster{1,t}.highcorr_numROIs_M > (0.25*length(tadcluster{1,t}.highcorr_numROIs_M)))

% find the intersection - ROIs that are highly correlated in all 3
% conditions
commonROIs = intersect(intersect(high_MS,high_V), high_M)

% plot the df/f0 of all ROIs highly correlated in all 3 conditions
figure;
plot(tadcluster{1,t}.alldff0(commonROIs,:)')
for i = 1:length(commonROIs)
    legend_commonROIs{i} = num2str(commonROIs(i))
end

legend(legend_commonROIs, 'Orientation', 'horizontal')
saveas(gcf, 'highly correlated ROIs all df_f0', 'png') 

%% find low correlation cells (1(self):10% of ROIs)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% copied from correlations_v2_xcorr lines 371-426
% based on multi df/f0
for t=1:length(tadcluster)
    if sum(size(tadcluster{1,t}.dff0_bystimtype{1,1})) > 0
        if length(tadcluster{1,t}.resp_ROIs) > 0
            for r = 1:size(tadcluster{1,t}.respROIdff0_maxR_sq_MS, 1)
                for c = 1:size(tadcluster{1,t}.respROIdff0_maxR_sq_MS, 2)
                    if tadcluster{1,t}.respROIdff0_maxR_sq_MS(r,c) < 0.2
                        highcorr(r, c) = 1;
                    else
                        highcorr(r, c) = 0;
                    end
                end
            end
            tadcluster{1,t}.respROIdff0_lowcorr_MS = logical(highcorr)
            clear('highcorr')
        end
    end
end            
            
% based on vis df/f0
for t=1:length(tadcluster)
    if sum(size(tadcluster{1,t}.dff0_bystimtype{1,1})) > 0
        if length(tadcluster{1,t}.resp_ROIs) > 0
            for r = 1:size(tadcluster{1,t}.respROIdff0_maxR_sq_V, 1)
                for c = 1:size(tadcluster{1,t}.respROIdff0_maxR_sq_V, 2)
                    if tadcluster{1,t}.respROIdff0_maxR_sq_V(r,c) < 0.2
                        highcorr(r, c) = 1;
                    else
                        highcorr(r, c) = 0;
                    end
                end
            end
            tadcluster{1,t}.respROIdff0_lowcorr_V = logical(highcorr)
            clear('highcorr')
        end
    end
end 
                
% based on mech df/f0
for t=1:length(tadcluster)
    if sum(size(tadcluster{1,t}.dff0_bystimtype{1,1})) > 0
        if length(tadcluster{1,t}.resp_ROIs) > 0
            for r = 1:size(tadcluster{1,t}.respROIdff0_maxR_sq_M, 1)
                for c = 1:size(tadcluster{1,t}.respROIdff0_maxR_sq_M, 2)
                    if tadcluster{1,t}.respROIdff0_maxR_sq_M(r,c) < 0.2
                        highcorr(r, c) = 1;
                    else
                        highcorr(r, c) = 0;
                    end
                end
            end
            tadcluster{1,t}.respROIdff0_lowcorr_M = logical(highcorr)
            clear('highcorr')
        end
    end
end 

% lines 594-604
for t = 1:length(tadcluster)
   if sum(size(tadcluster{1,t}.dff0_bystimtype{1,1})) > 0
        if length(tadcluster{1,t}.resp_ROIs) > 0
            %tadcluster{1,t}.lowcorr_numROIs_all = sum(tadcluster{1,t}.respROIdff0_highcorr);
            tadcluster{1,t}.lowcorr_numROIs_MS = sum(tadcluster{1,t}.respROIdff0_lowcorr_MS);
            tadcluster{1,t}.lowcorr_numROIs_V = sum(tadcluster{1,t}.respROIdff0_lowcorr_V);
            tadcluster{1,t}.lowcorr_numROIs_M = sum(tadcluster{1,t}.respROIdff0_lowcorr_M);
            %tadcluster{1,t}.lowcorr_numROIs_N = sum(tadcluster{1,t}.respROIdff0_highcorr_N);
        end
   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t=12
% find ROIs in tad 12 with at least 0.25*number of ROIs low correlated
low_MS = find(tadcluster{1,t}.lowcorr_numROIs_MS > (0.25*length(tadcluster{1,t}.highcorr_numROIs_MS)))
low_V = find(tadcluster{1,t}.lowcorr_numROIs_V > (0.25*length(tadcluster{1,t}.highcorr_numROIs_V)))
low_M = find(tadcluster{1,t}.lowcorr_numROIs_M > (0.25*length(tadcluster{1,t}.highcorr_numROIs_M)))

% find the intersection - ROIs that are highly correlated in all 3
% conditions
commonROIs_low = intersect(intersect(low_MS,low_V), low_M)

% plot the df/f0 of all ROIs low correlated in all 3 conditions
figure;
plot(tadcluster{1,t}.alldff0(commonROIs_low,:)')
for i = 1:length(commonROIs_low)
    legend_commonROIs_low{i} = num2str(commonROIs_low(i))
end

legend(legend_commonROIs_low, 'Orientation', 'horizontal')
saveas(gcf, 'low correlated ROIs all df_f0', 'png') 

%% Ok, but now we have ROIs with artifacts in them - let's remove that to get only ROIS that are nice looking
% remove any ROIs with > 2 df/f0
% highly correlated 
high_corr_ROIs = [];
for i = 1:length(commonROIs)
    if sum( find(tadcluster{1,t}.alldff0(commonROIs(i),:) > 1)) < 5
        high_corr_ROIs = [high_corr_ROIs, commonROIs(i)]
    end
end

% low correlated 
low_corr_ROIs = [];
for i = 1:length(commonROIs_low)
    if sum( find(tadcluster{1,t}.alldff0(commonROIs_low(i),:) > 1)) < 5
        low_corr_ROIs = [low_corr_ROIs, commonROIs_low(i)]
    end
end

% replot just the "good" ROIs from each
% plot the df/f0 of all ROIs highly correlated in all 3 conditions
figure;
plot(tadcluster{1,t}.alldff0(high_corr_ROIs,:)')
for i = 1:length(high_corr_ROIs)
    legend_high_corr_ROIs{i} = num2str(high_corr_ROIs(i))
end
legend(legend_high_corr_ROIs, 'Orientation', 'horizontal')
saveas(gcf, 'highly correlated ROIs all df_f0 good only', 'png')

figure;
plot(tadcluster{1,t}.alldff0(low_corr_ROIs,:)')
for i = 1:length(low_corr_ROIs)
    legend_low_corr_ROIs{i} = num2str(low_corr_ROIs(i))
end
legend(legend_low_corr_ROIs, 'Orientation', 'horizontal')
saveas(gcf, 'low correlated ROIs all df_f0 good only', 'png') 


%% Ok, now that I've done this for tad 12, do it for all tads so I can make histograms and do stat tests
for t = 1:length(tadcluster)
    if isfield(tadcluster{1,t}, 'highcorr_numROIs_MS')
    % find ROIs with at least 0.25*number of ROIs highly correlated
    high_MS = find(tadcluster{1,t}.highcorr_numROIs_MS > (0.25*length(tadcluster{1,t}.highcorr_numROIs_MS)))
    high_V = find(tadcluster{1,t}.highcorr_numROIs_V > (0.25*length(tadcluster{1,t}.highcorr_numROIs_V)))
    high_M = find(tadcluster{1,t}.highcorr_numROIs_M > (0.25*length(tadcluster{1,t}.highcorr_numROIs_M)))
    % find the intersection - ROIs that are highly correlated in all 3
    % conditions
    commonROIs = intersect(intersect(high_MS,high_V), high_M)
    
    % find ROIs in tad 12 with at least 0.25*number of ROIs low correlated
    low_MS = find(tadcluster{1,t}.lowcorr_numROIs_MS > (0.25*length(tadcluster{1,t}.highcorr_numROIs_MS)))
    low_V = find(tadcluster{1,t}.lowcorr_numROIs_V > (0.25*length(tadcluster{1,t}.highcorr_numROIs_V)))
    low_M = find(tadcluster{1,t}.lowcorr_numROIs_M > (0.25*length(tadcluster{1,t}.highcorr_numROIs_M)))
    % find the intersection - ROIs that are highly correlated in all 3
    % conditions
    commonROIs_low = intersect(intersect(low_MS,low_V), low_M)

    % remove any ROIs with > 2 df/f0
    % highly correlated 
    high_corr_ROIs = [];
    for i = 1:length(commonROIs)
        if sum( find(tadcluster{1,t}.alldff0(commonROIs(i),:) > 1)) < 5
            high_corr_ROIs = [high_corr_ROIs, commonROIs(i)]
        end
    end
    tadcluster{1,t}.high_corr_ROIs = high_corr_ROIs
    % low correlated 
    low_corr_ROIs = [];
    for i = 1:length(commonROIs_low)
        if sum( find(tadcluster{1,t}.alldff0(commonROIs_low(i),:) > 1)) < 5
            low_corr_ROIs = [low_corr_ROIs, commonROIs_low(i)]
        end
    end
    tadcluster{1,t}.low_corr_ROIs = low_corr_ROIs
    clear('high_MS', 'high_V', 'high_M', 'low_MS', 'low_V', 'low_M', 'commonROIs', 'commonROIs_low')
    end
end

%% Tectum shaped scatterplot marking high corr, low corr, other respROI and all ROIs
t = 12
for t = 1:length(tadcluster)
    if isfield(tadcluster{1,t}, 'highcorr_numROIs_MS')
        figure;
        hold on
        scatter(tadcluster{1,t}.ROIcenters(:,1), tadcluster{1,t}.ROIcenters(:,2), '+', 'MarkerEdgeColor', [0.5 0.5 0.5])
        scatter(tadcluster{1,t}.ROIcenters(tadcluster{1,t}.resp_ROIs,1), tadcluster{1,t}.ROIcenters(tadcluster{1,t}.resp_ROIs,2), 50, [0.7 0.7 0.7], 'filled')
        scatter(tadcluster{1,t}.ROIcenters(tadcluster{1,t}.high_corr_ROIs,1), tadcluster{1,t}.ROIcenters(tadcluster{1,t}.high_corr_ROIs, 2), 50, 'm', 'filled')
        scatter(tadcluster{1,t}.ROIcenters(tadcluster{1,t}.low_corr_ROIs,1), tadcluster{1,t}.ROIcenters(tadcluster{1,t}.low_corr_ROIs, 2), 50, 'g', 'filled')
        hold off
        annotation('textbox', 'Position', [0.14 0.21 .1 .1], 'String', ['All recorded ROIs (+)'], 'Color', [0.5 0.5 0.5], 'LineStyle', 'none' );
        annotation('textbox', 'Position', [0.14 0.16 .1 .1], 'String', ['All responding ROIs (o)'], 'Color', [0.7 0.7 0.7], 'LineStyle', 'none' );
        annotation('textbox', 'Position', [0.14 0.11 .1 .1], 'String', ['High correlation ROIs'], 'Color', 'm', 'LineStyle', 'none' );
        annotation('textbox', 'Position', [0.14 0.07 .1 .1], 'String', ['Low correlation ROIs'], 'Color', 'g', 'LineStyle', 'none' );
        title(sprintf('tad %d high and low correlation ROIs', t))
        fig_filename = sprintf('tad %d high and low correlation ROIs tectum shaped', t)
        saveas(gcf, fig_filename, 'png')
        close;
    end
end


   



%% What's different about these highly correlated cells - use all tads. 
% how do they compare to very low correlation cells?

%% first get data using code from correlations_v2_corrcoef

%% peaks by modality - mean and stddev of each stim type
% avg_peak is in tadcluster{1,t}.peak_avg

for t = 1:length(tadcluster)
    tadcluster{1,t}.stimmask = get_stimmask(tadcluster{1,t}.stimorder)
end
for t = 1:length(tadcluster)
    tadcluster{1,t}.peak_stddev_bymod = std_by_stimtype(tadcluster{1,t}.peak_bytrial, tadcluster{1,t}.stimmask)
end

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

%% then get values for high corr and low corr cells into 1 matrix
high_corr_data = [];
low_corr_data = [];
for t = 1:length(tadcluster)
    clear('tmp_data')
    if isfield(tadcluster{1,t}, 'high_corr_ROIs')
    if length(tadcluster{1,t}.high_corr_ROIs) > 0
    % high correlation
    tmp_data(:,1) = t * ones(1, length(tadcluster{1,t}.high_corr_ROIs)) %which tad
    tmp_data(:,2) = tadcluster{1,t}.high_corr_ROIs %which ROIs
    tmp_data(:,3:6) = tadcluster{1,t}.onset_time_mean(1:4,tadcluster{1,t}.high_corr_ROIs)' %onset time mean
    tmp_data(:,7:10) = tadcluster{1,t}.onset_time_stdev(1:4,tadcluster{1,t}.high_corr_ROIs)' %onset time stdev
    tmp_data(:, 11:14) = tadcluster{1,t}.peak_avg(1:4,tadcluster{1,t}.high_corr_ROIs)' %avg peak
    tmp_data(:, 15:18) = tadcluster{1,t}.peak_stddev_bymod(1:4,tadcluster{1,t}.high_corr_ROIs)' %avg peak
    tmp_data(:, 19) = tadcluster{1,t}.MSenh_peak(:,tadcluster{1,t}.high_corr_ROIs)'
    tmp_data(:, 20) = tadcluster{1,t}.MSEnh_numresponses(:,tadcluster{1,t}.high_corr_ROIs)'
    tmp_data(:, 21) = tadcluster{1,t}.unibias_numresponses(:,tadcluster{1,t}.high_corr_ROIs)'
    high_corr_data = [high_corr_data; tmp_data];
    clear('tmp_data')
    end
    if length(tadcluster{1,t}.low_corr_ROIs) > 0
    %low correlation
    tmp_data(:,1) = t * ones(1, length(tadcluster{1,t}.low_corr_ROIs)) %which tad
    tmp_data(:,2) = tadcluster{1,t}.low_corr_ROIs %which ROIs
    tmp_data(:,3:6) = tadcluster{1,t}.onset_time_mean(1:4,tadcluster{1,t}.low_corr_ROIs)' %onset time mean
    tmp_data(:,7:10) = tadcluster{1,t}.onset_time_stdev(1:4,tadcluster{1,t}.low_corr_ROIs)' %onset time stdev
    tmp_data(:, 11:14) = tadcluster{1,t}.peak_avg(1:4,tadcluster{1,t}.low_corr_ROIs)' %avg peak
    tmp_data(:, 15:18) = tadcluster{1,t}.peak_stddev_bymod(1:4,tadcluster{1,t}.low_corr_ROIs)' %avg peak
    tmp_data(:, 19) = tadcluster{1,t}.MSenh_peak(:,tadcluster{1,t}.low_corr_ROIs)'
    tmp_data(:, 20) = tadcluster{1,t}.MSEnh_numresponses(:,tadcluster{1,t}.low_corr_ROIs)'
    tmp_data(:, 21) = tadcluster{1,t}.unibias_numresponses(:,tadcluster{1,t}.low_corr_ROIs)'
    low_corr_data = [low_corr_data; tmp_data];
    clear('tmp_data')
    end
    end
end

%% Using the high and low corr data, compare each parameter

for i = 3:size(high_corr_data,2) %over all params
    [KS_H(i), KS_P(i)] = kstest2(high_corr_data(:,i), low_corr_data(:,i))
end

% what is significant?
sig_KS = find(KS_H == 1)

% 7 = onset time stdev MS

%% Compare high corr to all other responding ROIs

%% 1. Topographical distance

ROI_dist_high_corr = [];
ROI_dist_not_high_corr = [];
for t = 1:length(tadcluster)
    if isfield(tadcluster{1,t}, 'high_corr_ROIs')
        if length(tadcluster{1,t}.high_corr_ROIs) > 0
            other_respROIs = setdiff(tadcluster{1,t}.resp_ROIs, tadcluster{1,t}.high_corr_ROIs)
            for i = 1:length(other_respROIs)
                for j = 1:length(other_respROIs)
                    tmp_dist(i,j) = tadcluster{1,t}.ROIdist(other_respROIs(i),other_respROIs(j));
                end
            end
            ROI_dist_not_high_corr = [ROI_dist_not_high_corr; reshape(tmp_dist, [], 1)]
            clear('tmp_dist')
            for i = 1:length(tadcluster{1,t}.high_corr_ROIs)
                for j = 1:length(tadcluster{1,t}.high_corr_ROIs)
                    tmp_dist(i,j) = tadcluster{1,t}.ROIdist(tadcluster{1,t}.high_corr_ROIs(i),tadcluster{1,t}.high_corr_ROIs(j));
                end
            end
            ROI_dist_high_corr = [ROI_dist_high_corr; reshape(tmp_dist, [], 1)]
        end
    end
end
figure;
hist(ROI_dist_high_corr,20)
figure;
hist(ROI_dist_not_high_corr, 20)

[KS_H_topo, KS_P_topo] = kstest2(ROI_dist_high_corr, ROI_dist_not_high_corr)

figure;
h = histfit(ROI_dist_high_corr, 20, 'kernel')
figure;
h_others = histfit(ROI_dist_not_high_corr, 20, 'kernel')

figure;
hold on
h_others = histfit(ROI_dist_not_high_corr, 20, 'kernel')
h = histfit(ROI_dist_high_corr, 20, 'kernel')
hold off

figure;
%hold on
h_nothighcorr = hist(ROI_dist_not_high_corr, 20)
h_highcorr = hist(ROI_dist_high_corr, 20) 
hold off

%normalize to peak value
maxval_HC = max(h_highcorr)
norm_histHC = h_highcorr ./ maxval_HC
maxval_NHC = max(h_nothighcorr)
norm_histNHC = h_nothighcorr ./ maxval_NHC
figure;
bar(norm_histHC, 'histc', 'FaceColor', 'b')
figure;
bar(norm_histNHC, 'histc', 'FaceColor', 'g')

%% Ok, so this doesn't seem to capture the topography we see in "cluster" ROIs
% Investigate the ROIs that are in a cluster (cluster = ROI has a significant correlation with all other respROIs)
% they are stored in tadcluster{1,t}.correlated_ROIs_alldff0_common_AROI

% topographical distance within cluster vs out of all clusters
%eithin clusters
for t = 1:length(tadcluster)
    if isfield(tadcluster{1,t}, 'correlated_ROIs_alldff0_common_AROI')
    for i = 1:length(tadcluster{1,t}.correlated_ROIs_alldff0_common_AROI)
    if ~isempty(tadcluster{1,t}.correlated_ROIs_alldff0_common_AROI{i})
        for j = 1:length(tadcluster{1,t}.correlated_ROIs_alldff0_common_AROI{i})
            for k = 1:length(tadcluster{1,t}.correlated_ROIs_alldff0_common_AROI{i})
                tmp_data(j,k) = tadcluster{1,t}.ROIdist(tadcluster{1,t}.correlated_ROIs_alldff0_common_AROI{i}(j), tadcluster{1,t}.correlated_ROIs_alldff0_common_AROI{i}(k));
            end
        end
        tadcluster{1,t}.correlated_ROIs_alldff0_common_ROIdist{i} = tmp_data;
    end
    end
    end
end

% out of all clusters
for t = 1:length(tadcluster)
    if isfield(tadcluster{1,t}, 'correlated_ROIs_alldff0_common_AROI')
        in_rois = [];
        for i = 1:length(tadcluster{1,t}.correlated_ROIs_alldff0_common_AROI)
            in_rois = [in_rois; tadcluster{1,t}.correlated_ROIs_alldff0_common_AROI{i}(:)];
        end
        unique_in_rois = unique(in_rois);
        uncorr_rois = setdiff(tadcluster{1,t}.resp_ROIs, unique_in_rois);
        if length(uncorr_rois) > 0
        for m = 1:length(uncorr_rois)
            for n = 1:length(uncorr_rois)
                tmp_data(m,n) = tadcluster{1,t}.ROIdist(uncorr_rois(m), uncorr_rois(n));
            end
        end
        tadcluster{1,t}.uncorrelated_ROIs_alldff0_ROIdist{i} = tmp_data;
        clear('tmp_data')
        end
        clear('unique_in_rois', 'uncorr_rois')
    end
end

% combine all distances into 1 vector (and eliminate repeats from j,k and
% k,j)
% in cluster
roi_dist_correlatedROIs = [];
for t = 1:length(tadcluster)
    if isfield(tadcluster{1,t}, 'correlated_ROIs_alldff0_common_ROIdist')
        tmp_data = [];
        for i = 1:length(tadcluster{1,t}.correlated_ROIs_alldff0_common_ROIdist)
            if sum(size(tadcluster{1,t}.correlated_ROIs_alldff0_common_ROIdist{i})) > 1
                tmp_data = [tmp_data; reshape(tadcluster{1,t}.correlated_ROIs_alldff0_common_ROIdist{i}, [], 1)]
            end
        end
        roi_dist_correlatedROIs = [roi_dist_correlatedROIs; tmp_data];
    end
end
roi_dist_corrROIs_u = unique(roi_dist_correlatedROIs)

%out of cluster
roi_dist_uncorrelatedROIs = [];
for t = 1:length(tadcluster)
    if isfield(tadcluster{1,t}, 'uncorrelated_ROIs_alldff0_ROIdist')
        tmp_data = [];
        for i = 1:length(tadcluster{1,t}.uncorrelated_ROIs_alldff0_ROIdist)
            if sum(size(tadcluster{1,t}.uncorrelated_ROIs_alldff0_ROIdist{i})) > 1
                tmp_data = [tmp_data; reshape(tadcluster{1,t}.uncorrelated_ROIs_alldff0_ROIdist{i}, [], 1)];
            end
        end
        roi_dist_uncorrelatedROIs = [roi_dist_uncorrelatedROIs; tmp_data];
    end
end
roi_dist_uncorrROIs_u = unique(roi_dist_uncorrelatedROIs)

figure;
hist(roi_dist_uncorrROIs_u, 50)
figure;
hist(roi_dist_corrROIs_u, 50)

[KS_H_topo, KS_P_topo] = kstest2(roi_dist_uncorrROIs_u, roi_dist_corrROIs_u, 'Tail', 'smaller')
% significant, reject null hypothesis
figure;
hold on
ecdf(roi_dist_uncorrROIs_u)
ecdf(roi_dist_corrROIs_u)
hold off
legend('uncorr', 'corr')

%% Plot exp 12 tectum shaped scatterplot in a nicer way
% adapeted from correlations_v2_xcorr lines 540-578
t=12
rois = [];
list = [];
for i = 1:length(tadcluster{1,t}.correlated_ROIs_alldff0_common_AROI)
    rois = [rois tadcluster{1,t}.correlated_ROIs_alldff0_common_AROI{i}'] %all ROIs in all groups
    if ~isempty(tadcluster{1,t}.correlated_ROIs_alldff0_common_AROI{i})
        list = [list tadcluster{1,t}.resp_ROIs(i)]; %the ROIs we're using for the groups
    else
        continue
    end
end
xy_data = tadcluster{1,t}.ROIcenters(rois,:); %get ROI centers for all ROIs being included
% change size for each ROI-based cluster
sizes = [];
N = 1;
for j = 1:length(tadcluster{1,t}.correlated_ROIs_alldff0_common_AROI)
    if ~isempty(tadcluster{1,t}.correlated_ROIs_alldff0_common{j})
        len=length(tadcluster{1,t}.correlated_ROIs_alldff0_common{j});
        sizes = [sizes; N*ones(len,1)] ;
        N = N+1;
    end
end
plot_sizes = (1 ./ sizes) * 250;
%assign colors to each cluster
cmap = colormap(jet(N-1));
cluster_colors = cmap(sizes, :);
figure;
hold on
scatter(tadcluster{1,t}.ROIcenters(:,1), tadcluster{1,t}.ROIcenters(:,2), '+', 'MarkerEdgeColor', [0.7 0.7 0.7])
scatter(tadcluster{1,t}.ROIcenters(tadcluster{1,t}.resp_ROIs,1), tadcluster{1,t}.ROIcenters(tadcluster{1,t}.resp_ROIs,2), 50, 'filled', 'MarkerFaceColor', [0.5 0.5 0.5])
scatter(xy_data(:,1), xy_data(:,2), plot_sizes, cluster_colors, 'filled')
hold off
legend('unresponsive ROIs (+)', 'other responding ROIs', 'highly correlated ROIs')
