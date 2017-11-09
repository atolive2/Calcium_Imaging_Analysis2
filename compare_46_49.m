%% Compare stage 46 to 49

%%%%%%%% IMPORTANT NOTE: Half the data only has the smoothed df/f0
%%%%%%%% analysis, the other half has both raw and smoothed. Therefore,
%%%%%%%% when extracting information from allData, if the field "smoothed"
%%%%%%%% exists, then you need to use [param_name]_sm, not just
%%%%%%%% [param_name]

%% Import all the data
% start with xcorr analyzed data (not raw basics)

myFolder = 'D:/Torrey_calcium_imaging/compare_46-49/'; % May need to correct this.
if ~isdir(myFolder)
	errorMessage = sprintf('Error: The following folder does not exist:\n%s', myFolder);
	uiwait(warndlg(errorMessage));
	return;
end
filePattern = fullfile(myFolder, '*.mat');
matFiles = dir(filePattern)

counter = 0;
for k = 1:length(matFiles)
	matFilename = fullfile(myFolder, matFiles(k).name)
	matData = load(matFilename); % Retrieves a structure.
    
	% See if tadpole actually exists in the data structure.
	hasField = isfield(matData, 'tadcluster');
	if ~hasField
		% Alert user in popup window.
		warningMessage = sprintf('tadcluster is not in %s\n', matFilename);
		uiwait(warndlg(warningMessage));
		% It's not there, so skip the rest of the loop.
		continue; % Go to the next iteration.
    end
    for t = 1:length(matData.tadcluster)
        alltads{1,t+counter} = matData.tadcluster{1,t} % If you get to here, tadpole existed in the file.
        %counter = counter + 1
    end
    counter = length(alltads)
end

%% Now sort the tads by exp number, add stage and bathTBX info

% what tads do I have and in what order?
for t = 1:length(alltads)
    expnums(t) = alltads{1,t}.expnum
end
expnums_sort = sort(expnums)

% key - row1=expnum, row2=stage, row3=bathTBX 
key = [expnums_sort];

% sort data into new cell array by exp number, add stage and bathTBX fields
for t = 1:size(key,1)
    exp = key(t,1);
    col = find(expnums == exp)
    allData{1,t} = alltads{1,col};
    allData{1,t}.stage = key(t,2);
    allData{1,t}.bathRBX = key(t,3); %this is whether or not there is 0.1mM TBX (tubacurarine) while recording
end

% check to make sure that worked
for t = 1:length(alltads)
    expnums_AD(t) = allData{1,t}.expnum
end

% list the col num of all stage 46 and stage 49 tads separately
%this removes the one discobox 46
st46 = find(key(:,2) == 46);
st49 = find(key(:,2) == 49);
% generate x values for plotting purposes 
st46_X = ones(length(st46),1); % + (0.05 * rand(length(st46),1));
st49_X = 2*ones(length(st49),1);

%% Does the number of responding ROIs differ?

% gather count of respROIs and total ROI count into a matrix
for t = 1:length(allData)
    ct_respROIs(t,1) = length(allData{1,t}.resp_ROIs);
    ct_respROIs(t,2) = size(allData{1,t}.ROIcenters,1);
end

% scatterplot the values by stage - number of respROIs
figure;
hold on
scatter(st46_X, ct_respROIs(st46,1), 40, 'g', 'filled')
scatter(st49_X, ct_respROIs(st49,1), 40, [0.5 0 0.5], 'filled')
hold off
xlim([0.5 2.5])
ax = gca;
ax.XTick = [1 2]
ax.XTickLabel = [46 49]
title('Number of Responding ROIs by Exp')
xlabel('Stage')
ylabel('ROI count')
saveas(gcf, 'num resROIs by stage', 'png')

% scatterplot the values by stage - proportion of respROIs
figure;
hold on
scatter(st46_X, ct_respROIs(st46,2), 40, 'g', 'filled')
scatter(st49_X, ct_respROIs(st49,2), 40, [0.5 0 0.5], 'filled')
hold off
xlim([0.5 2.5])
ax = gca;
ax.XTick = [1 2];
ax.XTickLabel = [46 49];
title('Proportion of Responding ROIs by Exp')
xlabel('Stage')
ylabel('Proportion of ROIs that respond')
saveas(gcf, 'prop resROIs by stage', 'png')


%% Do 46s differ from 49s on any basic parameter?
% responses, area, peak, MSEnh, unimax_stimtype

%% calculate additional parameters
%% peaks by modality - mean and stddev of each stim type
% avg_peak is in allData{1,t}.peak_avg

for t = 1:length(allData)
    allData{1,t}.stimmask = get_stimmask(allData{1,t}.stimorder)
end
for t = 1:length(allData)
    if isfield(allData{1,t}, 'meanpeak_bytrial_sm')
        allData{1,t}.peak_stddev_bymod = std_by_stimtype(allData{1,t}.meanpeak_bytrial_sm, allData{1,t}.stimmask);
    else 
        allData{1,t}.peak_stddev_bymod = std_by_stimtype(allData{1,t}.peak_bytrial, allData{1,t}.stimmask);
    end 
end

%% jitter = variation in onset time (separate by modality)
% copy code from analyze_onset_time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%add stim onset time to allData%%%%
% exps 1-9 have onset at 0.5s 
% exps 10+ have onset at 2s
for t = 1:length(allData)
    if allData{1,t}.expnum <=9
        allData{1,t}.stim_onset = 0.5;
    elseif allData{1,t}.expnum >9
        allData{1,t}.stim_onset = 2;
    else
        fprintf('error exp %d, t')
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get time to 50% of the peak (= response onset time)
% same as used in analyze_onset_time
for t = 1:length(allData)
    if isfield(allData{1,t}, 'smoothed')
        for r = 1:size(allData{1,t}.smoothed,1) %over each ROI
            for tt = 1:size(allData{1,t}.smoothed,2) %over each trial
                half_peak = allData{1,t}.meanpeak_bytrial_sm(r, tt) / 2;
                tmp_onset = find((allData{1,t}.smoothed{r, tt} > half_peak) , 3);
                % get trial length (this takes care of exps with not 160 frames)
                trial_length = length(allData{1,t}.smoothed{r, tt});
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
    else
        for r = 1:size(allData{1,t}.df_f0,1) %over each ROI
            for tt = 1:size(allData{1,t}.df_f0,2) %over each trial
                half_peak = allData{1,t}.peak_bytrial(r, tt) / 2;
                tmp_onset = find((allData{1,t}.df_f0{r, tt} > half_peak) , 3);
                % get trial length (this takes care of exps with not 160 frames)
                trial_length = length(allData{1,t}.df_f0{r, tt});
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
    end
    
    % put onset_time into tadcluster
    allData{1,t}.onset_time = onset_time;    
end

% Get mean and sd of onset_time by ROI (this includes 0s - no response trials)
for t = 1:length(allData)
    allData{1,t}.onset_time_mean = mean_by_stimtype2(allData{1,t}.onset_time, allData{1,t}.stimmask);
    allData{1,t}.onset_time_stdev = std_by_stimtype(allData{1,t}.onset_time, allData{1,t}.stimmask);
end

%% multisensory enhancement (recalculate by num responses)

% use tadcluster{1,t}.MSEnh_peak 

% number of responses
for  t = 1:length(allData)
    if isfield(allData{1,t}, 'boolean_response_sm')
       for r = 1:size(allData{1,t}.boolean_response_sm, 1)
            multi = sum(allData{1,t}.boolean_response_sm(r, allData{1,t}.stimmask(:,1))) / sum(allData{1,t}.stimmask(:,1));
            vis = sum(allData{1,t}.boolean_response_sm(r, allData{1,t}.stimmask(:,2))) / sum(allData{1,t}.stimmask(:,2));
            mech = sum(allData{1,t}.boolean_response_sm(r, allData{1,t}.stimmask(:,3))) / sum(allData{1,t}.stimmask(:,3));
            if (vis + mech) == 0
                if multi > 0
                    allData{1,t}.MSEnh_numresponses(r) = 1;
                else
                    allData{1,t}.MSEnh_numresponses(r) = 0; 
                end
            else
                allData{1,t}.MSEnh_numresponses(r) = multi / (vis + mech);
            end
        end

    else
        
        for r = 1:size(allData{1,t}.boolean_response, 1)
            multi = sum(allData{1,t}.boolean_response(r, allData{1,t}.stimmask(:,1))) / sum(allData{1,t}.stimmask(:,1));
            vis = sum(allData{1,t}.boolean_response(r, allData{1,t}.stimmask(:,2))) / sum(allData{1,t}.stimmask(:,2));
            mech = sum(allData{1,t}.boolean_response(r, allData{1,t}.stimmask(:,3))) / sum(allData{1,t}.stimmask(:,3));
            if (vis + mech) == 0
                if multi > 0
                    allData{1,t}.MSEnh_numresponses(r) = 1;
                else
                    allData{1,t}.MSEnh_numresponses(r) = 0; 
                end
            else
                allData{1,t}.MSEnh_numresponses(r) = multi / (vis + mech);
            end
        end
    end
end

%% unisensory bias = ratio of vis responses to mech responses

% get number of responses to vis and mech, divide vis by mech
for  t = 1:length(allData)
    for r = 1:size(allData{1,t}.boolean_response, 1)
        vis = sum(allData{1,t}.boolean_response(r, allData{1,t}.stimmask(:,2))) / sum(allData{1,t}.stimmask(:,2));
        mech = sum(allData{1,t}.boolean_response(r, allData{1,t}.stimmask(:,3))) / sum(allData{1,t}.stimmask(:,3));
        if (vis + mech) == 0 %doesn't respond to either
            allData{1,t}.unibias_numresponses(r) = 0;
        else
            allData{1,t}.unibias_numresponses(r) = vis/(vis+mech);
        end
    end
end

%% Get unimax and stim type 
% using smoothed data, high/high (2/3) only

for t = 1:length(allData)
    if isfield(allData{1,t}, 'smoothed')
        [allData{1,t}.unimax_peakavg_sm, allData{1,t}.unimax_stimtype_sm] = max(allData{1,t}.peak_avg_sm(2:3,:));
    elseif ~isfield(allData{1,t}, 'unimax_peakavg')
        [allData{1,t}.unimax_peakavg, allData{1,t}.unimax_stimtype] = max(allData{1,t}.peak_avg(2:3,:));        
    end
end

%% combine data from responding ROIs from all exps into an array
idx = 1
for t = 1:length(allData)
    for r = 1:length(allData{1,t}.resp_ROIs)
        allRespROIs(idx, 1) = t;
        allRespROIs(idx, 2) = allData{1,t}.resp_ROIs(r);
        if isfield(allData{1,t}, 'smoothed')
            allRespROIs(idx,3:6) = allData{1,t}.area_avg_sm(1:4,allData{1,t}.resp_ROIs(r));
            allRespROIs(idx,7:10) = allData{1,t}.peak_avg_sm(1:4,allData{1,t}.resp_ROIs(r));
            allRespROIs(idx,11:14) = allData{1,t}.peakloc_avg_sm(1:4,allData{1,t}.resp_ROIs(r));
            allRespROIs(idx,15) = allData{1,t}.MSenh_peak_sm(:,allData{1,t}.resp_ROIs(r));
            allRespROIs(idx,16) = allData{1,t}.unimax_peakavg_sm(:,allData{1,t}.resp_ROIs(r));
            allRespROIs(idx,17) = allData{1,t}.unimax_stimtype_sm(:,allData{1,t}.resp_ROIs(r)); 
        else
            allRespROIs(idx,3:6) = allData{1,t}.area_avg(1:4,allData{1,t}.resp_ROIs(r));
            allRespROIs(idx,7:10) = allData{1,t}.peak_avg(1:4,allData{1,t}.resp_ROIs(r));
            allRespROIs(idx,11:14) = allData{1,t}.peakloc_avg(1:4,allData{1,t}.resp_ROIs(r));
            allRespROIs(idx,15) = allData{1,t}.MSenh_peak(:,allData{1,t}.resp_ROIs(r));
            allRespROIs(idx,16) = allData{1,t}.unimax_peakavg(:,allData{1,t}.resp_ROIs(r));
            allRespROIs(idx,17) = allData{1,t}.unimax_stimtype(:,allData{1,t}.resp_ROIs(r));
        
        end
        allRespROIs(idx, 18:21) = allData{1,t}.onset_time_mean(1:4, allData{1,t}.resp_ROIs(r));
        allRespROIs(idx, 22:25) = allData{1,t}.onset_time_stdev(1:4, allData{1,t}.resp_ROIs(r));
        allRespROIs(idx, 26) = allData{1,t}.MSEnh_numresponses(:, allData{1,t}.resp_ROIs(r));
        allRespROIs(idx, 27) = allData{1,t}.unibias_numresponses(:, allData{1,t}.resp_ROIs(r));
        allRespROIs(idx, 28) = allData{1,t}.stage;
        allRespROIs(idx, 29) = allData{1,t}.bathRBX;
        % num responses by stimtype
        % uni bias, onset time, MSEnh_numresponses
        % tad, roi, stage, bathTBX
       idx = idx +1; 
    end
end

%% Combining all cells, are there differences by stage in basic params?
st46 = find(allRespROIs(:,28) == 46);
st49 = find(allRespROIs(:,28) == 49);
labels = {'area MS'; 'area V'; 'area M'; 'area NS'; 'peak MS'; 'peak V'; 'peak M'; 'peak NS'; ... 
    'peak loc MS'; 'peak loc V'; 'peak loc M'; 'peak loc NS'; 'MSEnh peak'; 'unimax peak'; 'unimax stimtype'; ... 
    'onset time MS'; 'onset time V'; 'onset time M'; 'onset time NS';...
    'onset time SD MS'; 'onset time SD V'; 'onset time SD M'; 'onset time SD NS'; ...
     'MSEnh num resp'; 'Uni bias num resp'};
 
% make ECDF plots for each var 
for i = 3:27
    s46 = allRespROIs(st46, i);
    s49 = allRespROIs(st49, i);
    figure;
    hold on
    ecdf(s46)

    %set(h, 'Color', 'g')
    ecdf(s49)
    h = get(gca, 'children')
    set(h, 'LineWidth', 3)
    %j = get(gca, 'children')
    %set(j, 'LineWidth', 3)
    %set(j, 'Color', [0.5 0 0.5])
    hold off
    title(sprintf('%s ECDF by stage', labels{i-2}))
    xlabel(labels{i-2})
    ylabel('ROI count')
    fig_filename = sprintf('st 46 vs 49 ECDF of %s', labels{i-2})
    saveas(gcf, fig_filename, 'png')
    close;
end

% are any statistically different? (using kstest2)
for i = 3:27
    s46 = allRespROIs(st46, i);
    s49 = allRespROIs(st49, i);
    allRespROIs_H(i-2) = kstest2(s46, s49);
end
diff_vars = labels(find(allRespROIs_H))
%%%%%%%% sweet, almost everything is different


%% Does PCA differ by stage?
% copied from PCA_loadings
% This is based on Arseny's Elife paper plots
% looking for PCA loadings

% run PCA
% doc: https://www.mathworks.com/help/stats/pca.html

% get subset of allRespROIs data (eliminate area and peak loc)
allRespROIs_PCAsub = allRespROIs(:, [7:10, 15:27]);
labels_PCAsub = labels([5:8, 13:end])
% run PCA on all responding ROIs together
[coeff_all,score_all,latent_all,tsquared_all,explained_all,mu_all] = pca( [allRespROIs_PCAsub(:, :)]);
% run PCA on each stage separately
[coeff_46,score_46,latent_46,tsquared_46,explained_46,mu_46] = pca( [allRespROIs_PCAsub(st46, :)]);
[coeff_49,score_49,latent_49,tsquared_49,explained_49,mu_49] = pca( [allRespROIs_PCAsub(st49, :)]);

% the features (variables) put into allRespROIs(:, 3:27) are stored in labels

% plot coefficients and data onto PC1 x PC2 space
% https://www.mathworks.com/help/stats/biplot.html for documentation
figure;
biplot(coeff_all(:,1:2), 'Scores', score_all(:,1:2), 'VarLabels', labels_PCAsub)
title('PCA loadings all data')
figure;
biplot(coeff_46(:,1:2), 'Scores', score_46(:,1:2), 'VarLabels', labels_PCAsub)
title('PCA Loadings st 46 only')
figure;
biplot(coeff_49(:,1:2), 'Scores', score_49(:,1:2), 'VarLabels', labels_PCAsub)
title('PCA loadings st 49 only')

% plot amount of variance explained by each PC
figure;
bar(explained_all)
title('Variance explained by PC component all data')
xlabel('components')
ylabel('percent')
figure;
bar(explained_46)
title('Variance explained by PC component st46 only')
xlabel('components')
ylabel('percent')
figure;
bar(explained_49)
title('Variance explained by PC component st49 only')
xlabel('components')
ylabel('percent')

% plot MSEnh peak against PC1 value for all cells
figure;
scatter(score_all(:,1), allRespROIs_PCAsub(:,5))
title('MSEnh peak vs PC1 All Data')
xlabel('PC1')
ylabel('MSEnh peak')
saveas(gcf, 'MSEnh peak vs PC1 allRespROIs', 'png')

figure;
scatter(score_46(:,1), allRespROIs_PCAsub(st46,5))
title('MSEnh peak vs PC1 St46 only')
xlabel('PC1')
ylabel('MSEnh peak')
saveas(gcf, 'MSEnh peak vs PC1 st46 only', 'png')

figure;
scatter(score_49(:,1), allRespROIs_PCAsub(st49, 5))
title('MSEnh peak vs PC1 St49 only')
xlabel('PC1')
ylabel('MSEnh peak')
saveas(gcf, 'MSEnh peak vs PC1 st49 only', 'png')


%% Do 46s differ from 49s when assessing correlations?
% overall correlation, number of highcorr ROIs
% compare highcorr ROIs from each stim type
% this code is mostly copied from xcorr_assessment and applied to 46 and 49
% as found in allData. 

%% Define highcorr ROIs by actual ROI number for each modality (MS is done)


% First find significant overlap for MS, M and V
% significan overlap = at least 1/6*total num ROIS in correlated_ROIs_alldff0_int
clear('lens', 'int', 'roi_count')
for t = 1:length(allData)
    if sum(size(allData{1,t}.dff0_bystimtype{1,1})) > 0
        if length(allData{1,t}.resp_ROIs) > 0
            roi_count = (1/6)*length(allData{1,t}.resp_ROIs);
            for row = 1:size(allData{1,t}.correlated_ROIs_dff0_MS_int,1)
                for ct = 1:size(allData{1,t}.correlated_ROIs_dff0_MS_int,2)
                    lens(ct) = length(allData{1,t}.correlated_ROIs_dff0_MS_int{ct,row});
                end
                first_roi = find((lens > roi_count), 1)
                if isempty(first_roi) 
                    continue
                else
                    int = allData{1,t}.correlated_ROIs_dff0_MS_int{row, first_roi}

                    for col = first_roi:size(allData{1,t}.correlated_ROIs_dff0_MS_int,2)
                        if lens(col) > roi_count
                            int = intersect(int, allData{1,t}.correlated_ROIs_dff0_MS_int{row, col});
                        else
                            continue
                        end
                    end
                allData{1,t}.correlated_ROIs_dff0_MS_common{row} = int;
                end
                clear('lens', 'int')
            end
            clear('roi_count')
        end
    end
end


clear('lens', 'int', 'roi_count')
int = [];
for t = 1:length(allData)
    if sum(size(allData{1,t}.dff0_bystimtype{1,1})) > 0
        if length(allData{1,t}.resp_ROIs) > 0
            roi_count = (1/6)*length(allData{1,t}.resp_ROIs);
            for row = 1:size(allData{1,t}.correlated_ROIs_dff0_M_int,1)
                for ct = 1:size(allData{1,t}.correlated_ROIs_dff0_M_int,2)
                    lens(ct) = length(allData{1,t}.correlated_ROIs_dff0_M_int{ct,row});
                end
                first_roi = find((lens > roi_count), 1)
                
                if isempty(first_roi) 
                    continue
                else
                    int = allData{1,t}.correlated_ROIs_dff0_M_int{row, first_roi}

                    for col = first_roi:size(allData{1,t}.correlated_ROIs_dff0_M_int,2)
                        if lens(col) > roi_count
                            int = intersect(int, allData{1,t}.correlated_ROIs_dff0_M_int{row, col});
                        else
                            continue
                        end
                    end
                allData{1,t}.correlated_ROIs_dff0_M_common{row} = int;
                end
                clear('lens')
                int = [];
            end
            clear('roi_count')
        end
    end
end

clear('lens', 'int', 'roi_count')
for t = 1:length(allData)
    if sum(size(allData{1,t}.dff0_bystimtype{1,1})) > 0
        if length(allData{1,t}.resp_ROIs) > 0
            roi_count = (1/6)*length(allData{1,t}.resp_ROIs);
            for row = 1:size(allData{1,t}.correlated_ROIs_dff0_V_int,1)
                for ct = 1:size(allData{1,t}.correlated_ROIs_dff0_V_int,2)
                    lens(ct) = length(allData{1,t}.correlated_ROIs_dff0_V_int{ct,row});
                end
                first_roi = find((lens > roi_count), 1)
                if isempty(first_roi) 
                    continue
                elseif first_roi > length(allData{1,t}.resp_ROIs);
                    continue
                else
                    int = allData{1,t}.correlated_ROIs_dff0_V_int{row, first_roi}

                    for col = first_roi:size(allData{1,t}.correlated_ROIs_dff0_V_int,2)
                        if lens(col) > roi_count
                            int = intersect(int, allData{1,t}.correlated_ROIs_dff0_V_int{row, col});
                        else
                            continue
                        end
                    end
                allData{1,t}.correlated_ROIs_dff0_V_common{row} = int;
                end
                clear('lens', 'int')
            end
            clear('roi_count')
        end
    end
end

% Then take the ROI numbers and index them to the actual ROIs 
for t = 1:length(allData)
    if sum(size(allData{1,t}.dff0_bystimtype{1,1})) > 0
        if isfield(allData{1,t}, 'correlated_ROIs_dff0_MS_common')
            roi_list = allData{1,t}.resp_ROIs;
            for i = 1:length(allData{1,t}.correlated_ROIs_dff0_MS_common)
                if ~isempty(allData{1,t}.correlated_ROIs_dff0_MS_common(i))
                    allData{1,t}.correlated_ROIs_dff0_MS_common_AROI{i} = roi_list(allData{1,t}.correlated_ROIs_dff0_MS_common{i});
                end
            end
        end
        clear('roi_list')
    end
end

for t = 1:length(allData)
    if sum(size(allData{1,t}.dff0_bystimtype{1,1})) > 0
        if isfield(allData{1,t}, 'correlated_ROIs_dff0_M_common')
            roi_list = allData{1,t}.resp_ROIs;
            for i = 1:length(allData{1,t}.correlated_ROIs_dff0_M_common)
                if ~isempty(allData{1,t}.correlated_ROIs_dff0_M_common(i))
                    allData{1,t}.correlated_ROIs_dff0_M_common_AROI{i} = roi_list(allData{1,t}.correlated_ROIs_dff0_M_common{i});
                end
            end
        end
        clear('roi_list')
    end
end

for t = 1:length(allData)
    if sum(size(allData{1,t}.dff0_bystimtype{1,1})) > 0
        if isfield(allData{1,t}, 'correlated_ROIs_dff0_V_common')
            roi_list = allData{1,t}.resp_ROIs;
            for i = 1:length(allData{1,t}.correlated_ROIs_dff0_V_common)
                if ~isempty(allData{1,t}.correlated_ROIs_dff0_V_common(i))
                    allData{1,t}.correlated_ROIs_dff0_V_common_AROI{i} = roi_list(allData{1,t}.correlated_ROIs_dff0_V_common{i});
                end
            end
        end
        clear('roi_list')
    end
end


%% Find overlap in ROIs included in highcorr for each modality

% MS and M
for t = 1:length(allData)
    if isfield(allData{1,t}, 'correlated_ROIs_dff0_MS_common_AROI')
        if isfield(allData{1,t}, 'correlated_ROIs_dff0_M_common_AROI')
            M_rois = [];
            MS_rois = [];
            for i = 1:length(allData{1,t}.correlated_ROIs_dff0_M_common_AROI)
                M_rois = [M_rois; cell2mat(allData{1,t}.correlated_ROIs_dff0_M_common_AROI(i))];
            end
            for i = 1:length(allData{1,t}.correlated_ROIs_dff0_MS_common_AROI)
                MS_rois = [MS_rois; cell2mat(allData{1,t}.correlated_ROIs_dff0_MS_common_AROI(i))];
            end
            M_roisU = unique(M_rois)
            MS_roisU = unique(MS_rois)
            tmp = intersect(M_roisU, MS_roisU)
            common_ROIs_MS_M{1,t} =  tmp
            num_common_ROIs_MS_M(t, 1) = length(tmp);
            num_common_ROIs_MS_M(t, 2) = length(setdiff(MS_roisU, tmp));
            num_common_ROIs_MS_M(t, 3) = length(setdiff(M_roisU, tmp));
            num_common_ROIs_MS_M(t,4) = length(allData{1,t}.resp_ROIs);
            clear('tmp')
        else
           common_ROIs_MS_M{1,t} = [];
           num_common_ROIs_MS_M(t, :) = [NaN, NaN, NaN, length(allData{1,t}.resp_ROIs)]
        end
    else
        common_ROIs_MS_M{1,t} = [];
        num_common_ROIs_MS_M(t, :) = [NaN, NaN, NaN, length(allData{1,t}.resp_ROIs)]
    end
end

% MS and V
for t = 1:length(allData)
    if isfield(allData{1,t}, 'correlated_ROIs_dff0_MS_common_AROI')
        if isfield(allData{1,t}, 'correlated_ROIs_dff0_V_common_AROI')
            V_rois = [];
            MS_rois = [];
            for i = 1:length(allData{1,t}.correlated_ROIs_dff0_V_common_AROI)
                V_rois = [V_rois; cell2mat(allData{1,t}.correlated_ROIs_dff0_V_common_AROI(i))];
            end
            for i = 1:length(allData{1,t}.correlated_ROIs_dff0_MS_common_AROI)
                MS_rois = [MS_rois; cell2mat(allData{1,t}.correlated_ROIs_dff0_MS_common_AROI(i))];
            end
            V_roisU = unique(V_rois)
            MS_roisU = unique(MS_rois)
            tmp = intersect(V_roisU, MS_roisU)
            common_ROIs_MS_V{1,t} =  tmp
            num_common_ROIs_MS_V(t, 1) = length(tmp);
            num_common_ROIs_MS_V(t, 2) = length(setdiff(MS_roisU, tmp));
            num_common_ROIs_MS_V(t, 3) = length(setdiff(V_roisU, tmp));
            num_common_ROIs_MS_V(t,4) = length(allData{1,t}.resp_ROIs);
            clear('tmp')
        else
           common_ROIs_MS_V{1,t} = [];
           num_common_ROIs_MS_V(t, :) = [NaN, NaN, NaN, length(allData{1,t}.resp_ROIs)]
        end
    else
        common_ROIs_MS_V{1,t} = [];
        num_common_ROIs_MS_V(t, :) = [NaN, NaN, NaN, length(allData{1,t}.resp_ROIs)]
    end
end

% M and V
for t = 1:length(allData)
    if isfield(allData{1,t}, 'correlated_ROIs_dff0_M_common_AROI')
        if isfield(allData{1,t}, 'correlated_ROIs_dff0_V_common_AROI')
            V_rois = [];
            M_rois = [];
            for i = 1:length(allData{1,t}.correlated_ROIs_dff0_V_common_AROI)
                V_rois = [V_rois; cell2mat(allData{1,t}.correlated_ROIs_dff0_V_common_AROI(i))];
            end
            for i = 1:length(allData{1,t}.correlated_ROIs_dff0_M_common_AROI)
                M_rois = [M_rois; cell2mat(allData{1,t}.correlated_ROIs_dff0_M_common_AROI(i))];
            end
            V_roisU = unique(V_rois)
            M_roisU = unique(M_rois)
            tmp = intersect(V_roisU, M_roisU)
            common_ROIs_M_V{1,t} =  tmp
            num_common_ROIs_M_V(t, 1) = length(tmp);
            num_common_ROIs_M_V(t, 2) = length(setdiff(M_roisU, tmp));
            num_common_ROIs_M_V(t, 3) = length(setdiff(V_roisU, tmp));
            num_common_ROIs_M_V(t,4) = length(allData{1,t}.resp_ROIs);
            clear('tmp')
        else
           common_ROIs_M_V{1,t} = [];
           num_common_ROIs_M_V(t, :) = [NaN, NaN, NaN, length(allData{1,t}.resp_ROIs)]
        end
    else
        common_ROIs_M_V{1,t} = [];
        num_common_ROIs_M_V(t, :) = [NaN, NaN, NaN, length(allData{1,t}.resp_ROIs)]
    end
end

% id 46 vs 49
s46_tads = [];
s49_tads = [];
for t = 1:length(allData)
    if allData{1,t}.stage == 46
        s46_tads = [s46_tads, t];
    elseif  allData{1,t}.stage == 49
        s49_tads = [s49_tads, t];
    else
        continue
    end
end

% Plot stacked bar graph for each pairing for each stage
bar(num_common_ROIs_MS_M(s46_tads, 1:3), 'stacked')
legend({'both'; 'MS'; 'M'})
title('Stage 46 common Highcorr ROIs MS-M')
saveas(gcf, 'num commmon highcorr ROIs MS_M s46', 'png')

bar(num_common_ROIs_MS_M(s49_tads, 1:3), 'stacked')
legend({'both'; 'MS'; 'M'})
title('Stage 49 common Highcorr ROIs MS-M')
saveas(gcf, 'num commmon highcorr ROIs MS_M s49', 'png')

bar(num_common_ROIs_MS_V(s46_tads, 1:3), 'stacked')
legend({'both'; 'MS'; 'V'}, 'Location', 'northwest')
title('Stage 46 common Highcorr ROIs MS-V')
saveas(gcf, 'num commmon highcorr ROIs MS_V s46', 'png')

bar(num_common_ROIs_MS_V(s49_tads, 1:3), 'stacked')
legend({'both'; 'MS'; 'V'}, 'Location', 'northwest')
title('Stage 49 common Highcorr ROIs MS-V')
saveas(gcf, 'num commmon highcorr ROIs MS_V s49', 'png')

bar(num_common_ROIs_M_V(s46_tads, 1:3), 'stacked')
legend({'both'; 'M'; 'V'}, 'Location', 'northwest')
title('Stage 46 common Highcorr ROIs M-V')
saveas(gcf, 'num commmon highcorr ROIs M_V s46', 'png')

bar(num_common_ROIs_M_V(s49_tads, 1:3), 'stacked')
legend({'both'; 'M'; 'V'}, 'Location', 'northwest')
title('Stage 49 common Highcorr ROIs M-V')
saveas(gcf, 'num commmon highcorr ROIs M_V s49', 'png')

%% Combine across tadpoles to compare maxR vals across stim

%multi
tmp_data46 = [];
tmp_data49 = [];
for t = 1:length(allData)
    if ismember(t, s46_tads)
        if isfield(allData{1,t}, 'respROIdff0_maxR_sq_MS')
            for i = 1:size(allData{1,t}.respROIdff0_maxR_sq_MS,1)
                for j = 1:size(allData{1,t}.respROIdff0_maxR_sq_MS,2)
                    tmp_data46 = [tmp_data46, allData{1,t}.respROIdff0_maxR_sq_MS(i, (j+1):end)];
                end
            end                     
        end
    elseif ismember(t, s49_tads)
        if isfield(allData{1,t}, 'respROIdff0_maxR_sq_MS')
            for i = 1:size(allData{1,t}.respROIdff0_maxR_sq_MS,1)
                for j = 1:size(allData{1,t}.respROIdff0_maxR_sq_MS,2)
                    tmp_data49 = [tmp_data49, allData{1,t}.respROIdff0_maxR_sq_MS(i, (j+1):end)];
                end
            end                     
        end
    end
end
maxR_allrespROI_MS46 = tmp_data46;
maxR_allrespROI_MS49 = tmp_data49;

%vis
tmp_data46 = [];
tmp_data49 = [];
for t = 1:length(allData)
    if ismember(t, s46_tads)
        if isfield(allData{1,t}, 'respROIdff0_maxR_sq_V')
            for i = 1:size(allData{1,t}.respROIdff0_maxR_sq_V,1)
                for j = 1:size(allData{1,t}.respROIdff0_maxR_sq_V,2)
                    tmp_data46 = [tmp_data46, allData{1,t}.respROIdff0_maxR_sq_V(i, (j+1):end)];
                end
            end                     
        end
    elseif ismember(t, s49_tads)
        if isfield(allData{1,t}, 'respROIdff0_maxR_sq_V')
            for i = 1:size(allData{1,t}.respROIdff0_maxR_sq_V,1)
                for j = 1:size(allData{1,t}.respROIdff0_maxR_sq_V,2)
                    tmp_data49 = [tmp_data49, allData{1,t}.respROIdff0_maxR_sq_V(i, (j+1):end)];
                end
            end                     
        end
    end
end
maxR_allrespROI_V46 = tmp_data46;
maxR_allrespROI_V49 = tmp_data49;

%mech
tmp_data46 = [];
tmp_data49 = [];
for t = 1:length(allData)
    if ismember(t, s46_tads)
        if isfield(allData{1,t}, 'respROIdff0_maxR_sq_M')
            for i = 1:size(allData{1,t}.respROIdff0_maxR_sq_M,1)
                for j = 1:size(allData{1,t}.respROIdff0_maxR_sq_M,2)
                    tmp_data46 = [tmp_data46, allData{1,t}.respROIdff0_maxR_sq_M(i, (j+1):end)];
                end
            end                     
        end
    elseif ismember(t, s49_tads)
        if isfield(allData{1,t}, 'respROIdff0_maxR_sq_M')
            for i = 1:size(allData{1,t}.respROIdff0_maxR_sq_M,1)
                for j = 1:size(allData{1,t}.respROIdff0_maxR_sq_M,2)
                    tmp_data49 = [tmp_data49, allData{1,t}.respROIdff0_maxR_sq_M(i, (j+1):end)];
                end
            end                     
        end
    end
end
maxR_allrespROI_M46 = tmp_data46;
maxR_allrespROI_M49 = tmp_data49;

% make ECDF of all combined data
figure;
hold on
ecdf(maxR_allrespROI_MS46)
ecdf(maxR_allrespROI_V46)
ecdf(maxR_allrespROI_M46)
%ecdf(maxR_allrespROI_N)
title('ECDF of maxR st46 tads all respROIs')
xlabel('maxR')
ylabel('ROI proportion')
annotation('textbox', 'Position', [0.2 0.75 .1 .1], 'String', ['Multi'], 'Color', 'b', 'LineStyle', 'none' );
annotation('textbox', 'Position', [0.2 0.7 .1 .1], 'String', ['Vis'], 'Color', 'r', 'LineStyle', 'none' );
annotation('textbox', 'Position', [0.2 0.65 .1 .1], 'String', ['Mech'], 'Color', 'y', 'LineStyle', 'none' );
%annotation('textbox', 'Position', [0.2 0.6 .1 .1], 'String', ['No stim'], 'Color', [0.5 0 0.5], 'LineStyle', 'none' );
%annotation('textbox', 'Position', [0.2 0.55 .1 .1], 'String', 'N = 16, n = 474', 'Color', 'k', 'LineStyle', 'none');
hold off
saveas(gcf, 'ECDF of maxR st46 tads all respROIs', 'png')

figure;
hold on
ecdf(maxR_allrespROI_MS49)
ecdf(maxR_allrespROI_V49)
ecdf(maxR_allrespROI_M49)
%ecdf(maxR_allrespROI_N)
title('ECDF of maxR st49 tads all respROIs')
xlabel('maxR')
ylabel('ROI proportion')
annotation('textbox', 'Position', [0.2 0.75 .1 .1], 'String', ['Multi'], 'Color', 'b', 'LineStyle', 'none' );
annotation('textbox', 'Position', [0.2 0.7 .1 .1], 'String', ['Vis'], 'Color', 'r', 'LineStyle', 'none' );
annotation('textbox', 'Position', [0.2 0.65 .1 .1], 'String', ['Mech'], 'Color', 'y', 'LineStyle', 'none' );
%annotation('textbox', 'Position', [0.2 0.6 .1 .1], 'String', ['No stim'], 'Color', [0.5 0 0.5], 'LineStyle', 'none' );
%annotation('textbox', 'Position', [0.2 0.55 .1 .1], 'String', 'N = 16, n = 474', 'Color', 'k', 'LineStyle', 'none');
hold off
saveas(gcf, 'ECDF of maxR st49 tads all respROIs', 'png')

%% stat test all combined data
% Statistical difference across 46?
all_mods46 = [maxR_allrespROI_MS46', maxR_allrespROI_V46', maxR_allrespROI_M46'];
size(all_mods46)
[maxR_p46, maxR_tbl46, maxR_stats46] = kruskalwallis(all_mods46)
maxR_c46 = multcompare(maxR_stats46)
% MS diff from M and V, M and V not diff from each other

% Statistical difference across 49?
all_mods49 = [maxR_allrespROI_MS49', maxR_allrespROI_V49', maxR_allrespROI_M49'];
size(all_mods49)
[maxR_p49, maxR_tbl49, maxR_stats49] = kruskalwallis(all_mods49)
maxR_c49 = multcompare(maxR_stats49)
% MS, M and V all diff from each other

%% Highcorr vs not highcorr
% split by stage, MS-defined highcorr ROIs only

%% define all highcorr ROIs

tmp =  [];
for t = 1:length(allData)
    if isfield(allData{1,t}, 'correlated_ROIs_dff0_MS_common_AROI')
    for r = 1:length(allData{1,t}.correlated_ROIs_dff0_MS_common_AROI)
        tmp = [tmp; cell2mat(allData{1,t}.correlated_ROIs_dff0_MS_common_AROI(r))];
    end
    allData{1,t}.uniqueHighCorrROI =  unique(tmp)
    tmp = [];
    end
end


%% Are high corr ROIs more active overall?

highcorrROIs46 = [];
nothighcorrROIs46 = [];
highcorrROIs49 = [];
nothighcorrROIs49 = [];
for t = 1:length(allData)
    if ismember(t, s46_tads)
        if isfield(allData{1,t}, 'correlated_ROIs_dff0_MS_common_AROI')
            fprintf(num2str(t))
            for r = 1:length(allData{1,t}.resp_ROIs)
                if ismember(allData{1,t}.resp_ROIs(r), allData{1,t}.uniqueHighCorrROI)
                    highcorrROIs46 = [highcorrROIs46, allData{1,t}.sum_responses(allData{1,t}.resp_ROIs(r))];
                else
                    nothighcorrROIs46 = [nothighcorrROIs46, allData{1,t}.sum_responses(allData{1,t}.resp_ROIs(r))];
                end
            end
        end
    elseif ismember(t, s49_tads)
        if isfield(allData{1,t}, 'correlated_ROIs_dff0_MS_common_AROI')
            fprintf(num2str(t))
            for r = 1:length(allData{1,t}.resp_ROIs)
                if ismember(allData{1,t}.resp_ROIs(r), allData{1,t}.uniqueHighCorrROI)
                    highcorrROIs49 = [highcorrROIs49, allData{1,t}.sum_responses(allData{1,t}.resp_ROIs(r))];
                else
                    nothighcorrROIs49 = [nothighcorrROIs49, allData{1,t}.sum_responses(allData{1,t}.resp_ROIs(r))];
                end
            end
        end 
    else
        continue
    end
end

numHC46 = length(highcorrROIs46)
numNHC46 = length(nothighcorrROIs46)
numHC49 = length(highcorrROIs49)
numNHC49 = length(nothighcorrROIs49)

figure;
hold on 
ecdf(highcorrROIs46)
pause
ecdf(nothighcorrROIs46)
hold off
title('Total number of responses by high/not high corr st 46')

[h46, p46] = kstest2(highcorrROIs46, nothighcorrROIs46)

figure;
hold on 
ecdf(highcorrROIs49)
pause
ecdf(nothighcorrROIs49)
hold off
title('Total number of responses by high/not high corr st 49')

[h49, p49] = kstest2(highcorrROIs49, nothighcorrROIs49)

%% Is MSEnh_peak different between high corr and not high corr?

highcorrROIs46 = [];
nothighcorrROIs46 = [];
highcorrROIs49 = [];
nothighcorrROIs49 = [];
for t = 1:length(allData)
    if ismember(t, s46_tads)
        if isfield(allData{1,t}, 'correlated_ROIs_dff0_MS_common_AROI')
            fprintf(num2str(t))
            for r = 1:length(allData{1,t}.resp_ROIs)
                if ismember(allData{1,t}.resp_ROIs(r), allData{1,t}.uniqueHighCorrROI)
                    highcorrROIs46 = [highcorrROIs46, allData{1,t}.MSenh_peak(allData{1,t}.resp_ROIs(r))];
                else
                    nothighcorrROIs46 = [nothighcorrROIs46, allData{1,t}.MSenh_peak(allData{1,t}.resp_ROIs(r))];
                end
            end
        end
    elseif ismember(t, s49_tads)
        if isfield(allData{1,t}, 'correlated_ROIs_dff0_MS_common_AROI')
            fprintf(num2str(t))
            for r = 1:length(allData{1,t}.resp_ROIs)
                if ismember(allData{1,t}.resp_ROIs(r), allData{1,t}.uniqueHighCorrROI)
                    highcorrROIs49 = [highcorrROIs49, allData{1,t}.MSenh_peak(allData{1,t}.resp_ROIs(r))];
                else
                    nothighcorrROIs49 = [nothighcorrROIs49, allData{1,t}.MSenh_peak(allData{1,t}.resp_ROIs(r))];
                end
            end
        end 
    else
        continue
    end
end

numHC46 = length(highcorrROIs46)
numNHC46 = length(nothighcorrROIs46)
numHC49 = length(highcorrROIs49)
numNHC49 = length(nothighcorrROIs49)

figure;
hold on 
ecdf(highcorrROIs46)
pause
ecdf(nothighcorrROIs46)
hold off
title('MSEnh peak by high/not high corr st 46')

[h46, p46] = kstest2(highcorrROIs46, nothighcorrROIs46)

figure;
hold on 
ecdf(highcorrROIs49)
pause
ecdf(nothighcorrROIs49)
hold off
title('MSEnh peak by high/not high corr st 49')

[h49, p49] = kstest2(highcorrROIs49, nothighcorrROIs49)

%% Is average peak-multi different?

highcorrROIs46 = [];
nothighcorrROIs46 = [];
highcorrROIs49 = [];
nothighcorrROIs49 = [];
for t = 1:length(allData)
    if ismember(t, s46_tads)
        if isfield(allData{1,t}, 'correlated_ROIs_dff0_MS_common_AROI')
            fprintf(num2str(t))
            for r = 1:length(allData{1,t}.resp_ROIs)
                if ismember(allData{1,t}.resp_ROIs(r), allData{1,t}.uniqueHighCorrROI)
                    highcorrROIs46 = [highcorrROIs46, allData{1,t}.peak_avg(1, allData{1,t}.resp_ROIs(r))];
                else
                    nothighcorrROIs46 = [nothighcorrROIs46, allData{1,t}.peak_avg(1, allData{1,t}.resp_ROIs(r))];
                end
            end
        end
    elseif ismember(t, s49_tads)
        if isfield(allData{1,t}, 'correlated_ROIs_dff0_MS_common_AROI')
            fprintf(num2str(t))
            for r = 1:length(allData{1,t}.resp_ROIs)
                if ismember(allData{1,t}.resp_ROIs(r), allData{1,t}.uniqueHighCorrROI)
                    highcorrROIs49 = [highcorrROIs49, allData{1,t}.peak_avg(1, allData{1,t}.resp_ROIs(r))];
                else
                    nothighcorrROIs49 = [nothighcorrROIs49, allData{1,t}.peak_avg(1, allData{1,t}.resp_ROIs(r))];
                end
            end
        end 
    else
        continue
    end
end

numHC46 = length(highcorrROIs46)
numNHC46 = length(nothighcorrROIs46)
numHC49 = length(highcorrROIs49)
numNHC49 = length(nothighcorrROIs49)

figure;
hold on 
ecdf(highcorrROIs46)
pause
ecdf(nothighcorrROIs46)
hold off
title('Avg MS peak by high/not high corr st 46')

[h46, p46] = kstest2(highcorrROIs46, nothighcorrROIs46)

figure;
hold on 
ecdf(highcorrROIs49)
pause
ecdf(nothighcorrROIs49)
hold off
title('Avg MS peak by high/not high corr st 49')

[h49, p49] = kstest2(highcorrROIs49, nothighcorrROIs49)


%% Is average peak-unimax different?

highcorrROIs46 = [];
nothighcorrROIs46 = [];
highcorrROIs49 = [];
nothighcorrROIs49 = [];
for t = 1:length(allData)
    if ismember(t, s46_tads)
        if isfield(allData{1,t}, 'correlated_ROIs_dff0_MS_common_AROI')
            if isfield(allData{1,t}, 'peak_avg_sm')
                for r = 1:length(allData{1,t}.resp_ROIs)
                    if ismember(allData{1,t}.resp_ROIs(r), allData{1,t}.uniqueHighCorrROI)
                        highcorrROIs46 = [highcorrROIs46, max(allData{1,t}.peak_avg_sm(2:3, allData{1,t}.resp_ROIs(r)))];
                    else
                        nothighcorrROIs46 = [nothighcorrROIs46, max(allData{1,t}.peak_avg_sm(2:3, allData{1,t}.resp_ROIs(r)))];
                    end
                end
            else
                for r = 1:length(allData{1,t}.resp_ROIs)
                    if ismember(allData{1,t}.resp_ROIs(r), allData{1,t}.uniqueHighCorrROI)
                        highcorrROIs46 = [highcorrROIs46, max(allData{1,t}.peak_avg(2:3, allData{2:3,t}.resp_ROIs(r)))];
                    else
                        nothighcorrROIs46 = [nothighcorrROIs46, max(allData{1,t}.peak_avg(2:3, allData{1,t}.resp_ROIs(r)))];
                    end
                end
            end
        end
    elseif ismember(t, s49_tads)
        if isfield(allData{1,t}, 'correlated_ROIs_dff0_MS_common_AROI')
            if isfield(allData{1,t}, 'peak_avg_sm')
                for r = 1:length(allData{1,t}.resp_ROIs)
                    if ismember(allData{1,t}.resp_ROIs(r), allData{1,t}.uniqueHighCorrROI)
                        highcorrROIs49 = [highcorrROIs49, max(allData{1,t}.peak_avg_sm(2:3, allData{1,t}.resp_ROIs(r)))];
                    else
                        nothighcorrROIs49 = [nothighcorrROIs49, max(allData{1,t}.peak_avg_sm(2:3, allData{1,t}.resp_ROIs(r)))];
                    end
                end
            else
                for r = 1:length(allData{1,t}.resp_ROIs)
                    if ismember(allData{1,t}.resp_ROIs(r), allData{1,t}.uniqueHighCorrROI)
                        highcorrROIs49 = [highcorrROIs49, max(allData{1,t}.peak_avg(2:3, allData{1,t}.resp_ROIs(r)))];
                    else
                        nothighcorrROIs49 = [nothighcorrROIs49, max(allData{1,t}.peak_avg(2:3, allData{1,t}.resp_ROIs(r)))];
                    end
                end 
            
            end
        end 
    else
        continue
    end
end

numHC46 = length(highcorrROIs46)
numNHC46 = length(nothighcorrROIs46)
numHC49 = length(highcorrROIs49)
numNHC49 = length(nothighcorrROIs49)

figure;
hold on 
ecdf(highcorrROIs46)
pause
ecdf(nothighcorrROIs46)
hold off
title('Avg unimax peak by high/not high corr st 46')

figure;
hold on 
ecdf(highcorrROIs49)
pause
ecdf(nothighcorrROIs49)
hold off
title('Avg unimax peak by high/not high corr st 49')

[h46, p46] = kstest2(highcorrROIs46, nothighcorrROIs46)
[h49, p49] = kstest2(highcorrROIs49, nothighcorrROIs49)

%% N and n values
how_many46 = [];
how_many49 = [];
for t = 1:length(allData)
    tmp = [t, length(allData{1,t}.resp_ROIs), length(allData{1,t}.sum_responses), length(allData{1,t}.stimorder)];
    if allData{1,t}.stage == 46
        how_many46 = [how_many46; tmp];
    elseif allData{1,t}.stage == 49
        how_many49 = [how_many49; tmp];
    else
        fprintf('%s', t)
        continue
    end
end


%% Overall variability
% why are preps different?

%% Visual responsiveness

% how responsive is a given prep?
for t = 1:length(allData)
    V_tr = find(allData{1,t}.stimorder == 2)
    if isfield(allData{1,t}, 'boolean_response_sm')
        allData{1,t}.resp_ctV =  sum(allData{1,t}.boolean_response_sm(:, V_tr), 2);
    else
        allData{1,t}.resp_ctV =  sum(allData{1,t}.boolean_response(:, V_tr), 2);
    end
    clear('V_tr')
end

% Plot proportion of responses over trials * ROIs for 46 vs 49
to_plot46 = [];
to_plot49 = [];
for t = 1:length(allData)
    total_trV = length(allData{1,t}.stimorder) * length(allData{1,t}.resp_ctV);
    if ismember(t, s46_tads)
        to_plot46 = [to_plot46, (sum(allData{1,t}.resp_ctV) / total_trV)];
    elseif ismember(t, s49_tads)
        to_plot49 = [to_plot49, (sum(allData{1,t}.resp_ctV) / total_trV)];
    end
end

figure;
hold on
scatter(st46_X, to_plot46)
scatter(st49_X, to_plot49)
scatter(1, median(to_plot46), 50, 'filled')
scatter(2, median(to_plot49), 50, 'filled')
hold off
xlim([0.5 2.5])
ax = gca;
ax.XTick = [1 2]
ax.XTickLabel = [46 49]
title('Proportion of trials*ROIs with V response by stage')
xlabel('Stage')
ylabel('Proportion responses')
saveas(gcf, 'Proportion of trials x ROIs with V response by stage', 'png')

[h, p] = kstest2(to_plot46, to_plot49)

%% Is visual responsiveness correlated with highcorr ROIs?
highcorrROIs46 = [];
nothighcorrROIs46 = [];
highcorrROIs49 = [];
nothighcorrROIs49 = [];
for t = 1:length(allData)
    if ismember(t, s46_tads)
        if isfield(allData{1,t}, 'correlated_ROIs_dff0_MS_common_AROI')
            for r = 1:length(allData{1,t}.resp_ROIs)
                if ismember(allData{1,t}.resp_ROIs(r), allData{1,t}.uniqueHighCorrROI)
                    highcorrROIs46 = [highcorrROIs46, allData{1,t}.resp_ctV(allData{1,t}.resp_ROIs(r))];
                else
                    nothighcorrROIs46 = [nothighcorrROIs46, allData{1,t}.resp_ctV(allData{1,t}.resp_ROIs(r))];
                end
            end
        end
    elseif ismember(t, s49_tads)
        if isfield(allData{1,t}, 'correlated_ROIs_dff0_MS_common_AROI')
            for r = 1:length(allData{1,t}.resp_ROIs)
                if ismember(allData{1,t}.resp_ROIs(r), allData{1,t}.uniqueHighCorrROI)
                    highcorrROIs49 = [highcorrROIs49, allData{1,t}.resp_ctV(allData{1,t}.resp_ROIs(r))];
                else
                    nothighcorrROIs49 = [nothighcorrROIs49, allData{1,t}.resp_ctV(allData{1,t}.resp_ROIs(r))];
                end
            end
        end 
    else
        continue
    end
end

figure;
hold on 
ecdf(highcorrROIs46)
pause
ecdf(nothighcorrROIs46)
hold off
title('Visual responsiveness by high/not high corr st 46')

[h46, p46] = kstest2(highcorrROIs46, nothighcorrROIs46)

figure;
hold on 
ecdf(highcorrROIs49)
pause
ecdf(nothighcorrROIs49)
hold off
title('Visual responsiveness by high/not high corr st 49')

[h49, p49] = kstest2(highcorrROIs49, nothighcorrROIs49)

%% How many highcorr ROIs per tadpole?
% as a proportion of the total number of responding ROIs

% Plot proportion of responses over trials * ROIs for 46 vs 49
to_plot46 = [];
to_plot49 = [];
for t = 1:length(allData)
    if ismember(t, s46_tads)
        if isfield( allData{1,t}, 'uniqueHighCorrROI')
            to_plot46 = [to_plot46, (length(allData{1,t}.uniqueHighCorrROI) / length(allData{1,t}.resp_ROIs))];
        else
            to_plot46 = [to_plot46, 0];
        end
    elseif ismember(t, s49_tads)
        if isfield( allData{1,t}, 'uniqueHighCorrROI')
            to_plot49 = [to_plot49, (length(allData{1,t}.uniqueHighCorrROI) / length(allData{1,t}.resp_ROIs))];
        else
            to_plot49 = [to_plot49, 0];
        end
    end
end

figure;
hold on
scatter(st46_X, to_plot46)
scatter(st49_X, to_plot49)
scatter(1, median(to_plot46), 50, 'filled')
scatter(2, median(to_plot49), 50, 'filled')
hold off
xlim([0.5 2.5])
ax = gca;
ax.XTick = [1 2]
ax.XTickLabel = [46 49]
title('Proportion of ROIs that are highcorr by stage')
xlabel('Stage')
ylabel('Proportion ROIs')
saveas(gcf, 'Proportion of ROIs that are highcorr by stage', 'png')

[h, p] = kstest2(to_plot46, to_plot49)











