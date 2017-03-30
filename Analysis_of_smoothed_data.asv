%% Analysis of "Good" Experiments
% This code will import all "basic analysis" files, then smooth the df/f0
% data and build on it. 
% Steps needed
    % import tadpole{1, :}
    % smooth df/f0
    % recalulate all extracted values on smoothed data
    % eliminate any poor/annoying exps (e.g. ones without 160
        % frames/trial or 500x500 image size)
    % determine "responding" cells for further analysis
    % calculate difference in peak latency for multisensory vs unisensory
    % calculate average peak for each category

%% Import tadpole struct for all experiments
myFolder = 'F:/Calcium_Imaging_Analysis/analyzed_compiled/Smoothed_analysis/'; % May need to correct this.
%mkdir([myFolder 'figures']);
if ~isdir(myFolder)
	errorMessage = sprintf('Error: The following folder does not exist:\n%s', myFolder);
	uiwait(warndlg(errorMessage));
	return;
end
filePattern = fullfile(myFolder, 'exp*.mat');
matFiles = dir(filePattern)

for k = 1:length(matFiles)
	matFilename = fullfile(myFolder, matFiles(k).name)
	matData = load(matFilename); % Retrieves a structure.
    
	% See if tadpole actually exists in the data structure.
	hasField = isfield(matData, 'tadpole');
	if ~hasField
		% Alert user in popup window.
		warningMessage = sprintf('tadpole is not in %s\n', matFilename);
		uiwait(warndlg(warningMessage));
		% It's not there, so skip the rest of the loop.
		continue; % Go to the next iteration.
	end
	tadpole{1,k} = matData.tadpole % If you get to here, tadpole existed in the file.
    tadpole{1,k}.matFilename = matFilename
end

%% Smooth df/f0

for t = 1:length(tadpole)
    for i = 1:size(tadpole{1,t}.df_f0, 1)
        for j = 1:size(tadpole{1,t}.df_f0, 2)
            tadpole{1,t}.smoothed{i,j} = smooth(tadpole{1,t}.df_f0{i,j}(:,:), 8, 'moving');
        end
    end
end

%% Recalculate all extracted values using smoothed data

for t = 1:length(tadpole)
    % Extract parameters for each trial
    [ tadpole{1,t}.area_bytrial_sm ] = calc_area( tadpole{1,t}.smoothed, 46 )
    [ tadpole{1,t}.meanpeak_bytrial_sm, tadpole{1,t}.peakloc_bytrial_sm, tadpole{1,t}.peak_bytrial_sm ] = calc_peak2( tadpole{1,t}.smoothed, 5, 5)

    % Define response/no response
    [ tadpole{1,t}.boolean_response_sm, tadpole{1,t}.sum_responses_sm ] = get_respondingROIs2( tadpole{1,t}.area_bytrial_sm, tadpole{1,t}.meanpeak_bytrial_sm, tadpole{1,t}.peakloc_bytrial_sm )

    % Find average area, peak and peakloc for each ROI for each stim type
    [ tadpole{1,t}.area_avg_sm ] = mean_by_stimtype ( tadpole{1,t}.area_bytrial_sm, tadpole{1,t}.stimmask )
    [ tadpole{1,t}.peak_avg_sm ] = mean_by_stimtype2 ( tadpole{1,t}.meanpeak_bytrial_sm, tadpole{1,t}.stimmask )
    [ tadpole{1,t}.peakloc_avg_sm ] = mean_by_stimtype2 ( tadpole{1,t}.peakloc_bytrial_sm, tadpole{1,t}.stimmask )

    % calculate MS enhancement
    % this is only on high/high presentations (1,2,3) because argh
    % calc_MSehancement
    [ tadpole{1,t}.MSenh_area_sm ] = calc_MSenhancement( tadpole{1,t}.area_avg_sm )
    [ tadpole{1,t}.MSenh_peak_sm ] = calc_MSenhancement( tadpole{1,t}.peak_avg_sm )
    [ tadpole{1,t}.MSenh_peakloc_sm ] = calc_MSenhancement( tadpole{1,t}.peakloc_avg_sm )

    % collect all values for a single stim type into a single place for easy
    % histogram - ing.
    [ tadpole{1,t}.stim_vals_area_sm ] = get_all_bystim( tadpole{1,t}.stimmask, unique(tadpole{1,t}.stimorder), tadpole{1,t}.area_bytrial_sm )
    [ tadpole{1,t}.stim_vals_meanpeak_sm ] = get_all_bystim( tadpole{1,t}.stimmask, unique(tadpole{1,t}.stimorder), tadpole{1,t}.meanpeak_bytrial_sm )
    [ tadpole{1,t}.stim_vals_peakloc_sm ] = get_all_bystim( tadpole{1,t}.stimmask, unique(tadpole{1,t}.stimorder), tadpole{1,t}.peakloc_bytrial_sm )
end

%% Select responding ROIs

% What is the distribution of peak values?
peak_vals_toplot = [];
for t = 1:length(tadpole)
    peak_vals_toplot = [peak_vals_toplot; [tadpole{1,t}.meanpeak_bytrial_sm(:)]];
end
%peakvals_toplot = cell2mat(peak_vals_toplot);
hist(peak_vals_toplot, 1000)
%axes([0 100 0 10])

% find trials with peak bigger than 5 and eliminate
for t = 1:length(tadpole)
    tadpole{1,t}.goodTrials = ones(size(tadpole{1,t}.meanpeak_bytrial_sm));
    for i = 1:size(tadpole{1,t}.meanpeak_bytrial_sm,1)
        for j = 1:size(tadpole{1,t}.meanpeak_bytrial_sm,2)
            if tadpole{1,t}.meanpeak_bytrial_sm(i,j) > 5
                tadpole{1,t}.goodROIs(i,j) = 0;
            end
        end
    end
end

% Find Good ROIs (respond at least 4 times = peak > 0.15)
for t = 1:length(tadpole)
    %clear('tadpole{1,t}.goodROIs')
    tadpole{1,t}.goodROIs = ones(size(tadpole{1,t}.meanpeak_bytrial_sm, 1), 1);
    
    for i = 1:size(tadpole{1,t}.meanpeak_bytrial_sm,1)
        tmpdata = tadpole{1,t}.meanpeak_bytrial_sm(i,:);
        findrois = find(tmpdata > 0.1) % && (tmpdata < 5));
            if length(findrois) < 4
                tadpole{1,t}.goodROIs(i) = 0;
            end
        clear('tmpdata', 'findrois')
    end
end
% convert goodROIs to logical (for use later)
for t = 1:length(tadpole)
    tadpole{1,t}.goodROIslog = logical(tadpole{1,t}.goodROIs)
end

% how many good ROIs per experiment?
for t = 1:length(tadpole)
    num_goodROIs(t) = sum(tadpole{1,t}.goodROIs)
end
% what proportion of rois are good?
for t = 1:length(tadpole)
    prop_goodROIs(t) = num_goodROIs(t)/length(tadpole{1,t}.goodROIs)
end
% what experiments are good (e.g. proportion of ROIs that respond >
% 0.2/20%)
goodExp = prop_goodROIs > 0.2
% 9 good experiments - let's just use those. 

%% Analyse responding ROIs of good experiments only
% good experiments: if goodExp(t)
% good ROIs in good exps: if tadpole{1,t}.goodROIslog
idx = 1
for t = 1:length(tadpole)
    if goodExp(t)
        for r = 1:length(tadpole{1,t}.somaticROIs)
            if tadpole{1,t}.goodROIslog(r)
                goodtadpole(idx).exp = tadpole{1,t}.expnum;
                goodtadpole(idx).roi = r;
                goodtadpole(idx).smoothed_df_f0 = tadpole{1,t}.smoothed(r,:);
                goodtadpole(idx).area_bytrial = cell2mat(tadpole{1,t}.area_bytrial_sm(r, :));
                goodtadpole(idx).meanpeak_bytrial = tadpole{1,t}.meanpeak_bytrial_sm(r, :);
                goodtadpole(idx).peakloc_bytrial = tadpole{1,t}.peakloc_bytrial_sm(r, :);
                goodtadpole(idx).sum_responses = tadpole{1,t}.sum_responses_sm(r);
                goodtadpole(idx).area_avg = tadpole{1,t}.area_avg_sm(:, r);
                goodtadpole(idx).peak_avg = tadpole{1,t}.peak_avg_sm(:, r);
                goodtadpole(idx).peakloc_avg = tadpole{1,t}.peakloc_avg_sm(:, r);
                goodtadpole(idx).MSenh_peak = tadpole{1,t}.MSenh_peak_sm(r);
                goodtadpole(idx).MSenh_peakloc = tadpole{1,t}.MSenh_peakloc_sm(r);
                goodtadpole(idx).unimax_peakavg = tadpole{1,t}.unimax_peakavg_sm;
                goodtadpole(idx).unimax_stimtype = tadpole{1,t}.unimax_stimtype_sm;
                goodtadpole(idx).multimax_peakavg = tadpole{1,t}.multimax_peakavg_sm;
                goodtadpole(idx).stimorder = tadpole{1,t}.stimorder;
                goodtadpole(idx).roiloc = tadpole{1,t}.somaticROICenters{1,r}(1).Centroid;
                idx = idx + 1
            end
            
        end
    end
end


%% Analysis of peak variability
% calculate difference in peak latency for multisensory vs unisensory
% calculate average peak latency for each category

for r = 1:length(goodtadpole) % over all rois
    uniquestims = unique(goodtadpole(r).stimorder);
    for s = 1:length(uniquestims) % over each stim type
        stimidx = uniquestims(s);
        trials_touse = find(goodtadpole(r).stimorder == stimidx);
        avg_peakloc(r, s) = mean(goodtadpole(r).peakloc_bytrial(1, trials_touse));
        std_peakloc(r, s) = std(goodtadpole(r).peakloc_bytrial(1, trials_touse));
    end
end

figure;
errorbar(avg_peakloc(:,1), std_peakloc(:,1))
figure;
errorbar(avg_peakloc(:,2), std_peakloc(:,2))
figure;
errorbar(avg_peakloc(:,3), std_peakloc(:,3))

figure;
hold on
errorbar(avg_peakloc(:,1), std_peakloc(:,1))
errorbar(avg_peakloc(:,2), std_peakloc(:,2))
errorbar(avg_peakloc(:,3), std_peakloc(:,3))
errorbar(avg_peakloc(:,4), std_peakloc(:,4))
hold off
legend ('multi', 'vis', 'mech', 'no stim', 'orientation', 'horizontal')

% 1 way anova with multiple conparisons on multi/vis/mech/no stim
[ p, tbl, stats ] = anova1(avg_peakloc(:,1:4))
[results, means] = multcompare(stats,'CType','bonferroni')

% run 1 way anova on uni versus multi only
uniavg = mean(avg_peakloc(:,2:3)');
totestdata = [avg_peakloc(:,1) uniavg']
[ p_um, tbl_um, stats_um ] = anova1(totestdata)
% p = 0.06 so not significant. 

% figure
multi_X = 0.7 + rand(1,319)*(1.3-0.7)
vis_X = 1.7 + rand(1,319)*(2.3-1.7)
mech_X = 2.7 + rand(1,319)*(3.3-2.7)
none_X = 3.7 + rand(1,319)*(4.3-3.7)
% plot means with stdev of means
avgof_avgpeakloc = mean(avg_peakloc) 
stdof_avgpeakloc = std(avg_peakloc) 

figure;
hold on
% all means as circles
scatter(multi_X, avg_peakloc(:,1))
scatter(vis_X, avg_peakloc(:,2))
scatter(mech_X, avg_peakloc(:,3))
scatter(none_X, avg_peakloc(:,4))
errorbar(avgof_avgpeakloc(1:4), stdof_avgpeakloc(1:4), 'ob', 'linewidth', 2, 'Markersize', 10, 'MarkerFaceColor', 'b') % mean with stdev
hold off
title('Mean Peak Location')
ylabel('Peak Location in seconds')
ax=gca;
xsize = 160
ax.YTick = [0, xsize/7, (xsize/7)*2, (xsize/7)*3, (xsize/7)*4, (xsize/7)*5, (xsize/7)*6, (xsize/7)*7];
ax.YTickLabel = {'0','1', '2', '3', '4', '5', '6', '7'};
xlabel('Stimulus Type')
ax.XTick = [1, 2, 3, 4];
ax.XTickLabel = {'Multi','Visual', 'Mech', 'None'};
fig_filename='Mean Peak Location RespROIs only';
saveas(gcf,fig_filename,'png');

% calculate difference in variation within trials of a single ROI
