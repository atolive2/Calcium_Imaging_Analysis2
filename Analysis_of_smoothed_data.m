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
myFolder = 'D:\\Torrey_calcium_imaging\newData_20170924' %'F:/Calcium_Imaging_Analysis/analyzed_compiled/Smoothed_analysis/'; % May need to correct this.
%mkdir([myFolder 'figures']);
if ~isdir(myFolder)
	errorMessage = sprintf('Error: The following folder does not exist:\n%s', myFolder);
	uiwait(warndlg(errorMessage));
	return;
end
filePattern = fullfile(myFolder, 'tadpole*.mat');
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
    [ tadpole{1,t}.area_avg_sm ] = mean_by_stimtype2 ( cell2mat(tadpole{1,t}.area_bytrial_sm), tadpole{1,t}.stimmask )
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

% Is the standard deviation (variability within ROIs) different?

% 1 way anova with multiple conparisons on multi/vis/mech/no stim
[ p_std, tbl_std, stats_std ] = anova1(std_peakloc(:,1:4))
[results_std, means_std] = multcompare(stats_std,'CType','bonferroni')

% run 1 way anova on uni versus multi only
uniavg_std = mean(std_peakloc(:,2:3)');
totestdata_std = [std_peakloc(:,1) uniavg']
[ p_um_std, tbl_um_std, stats_um_std ] = anova1(totestdata_std)
% p = 0.06 so not significant. 

% figure
multi_X = 0.7 + rand(1,319)*(1.3-0.7)
vis_X = 1.7 + rand(1,319)*(2.3-1.7)
mech_X = 2.7 + rand(1,319)*(3.3-2.7)
none_X = 3.7 + rand(1,319)*(4.3-3.7)
% plot means with stdev of means
avgof_stdpeakloc = mean(std_peakloc) 
stdof_stdpeakloc = std(std_peakloc) 

figure;
hold on
% all means as circles
scatter(multi_X, std_peakloc(:,1))
scatter(vis_X, std_peakloc(:,2))
scatter(mech_X, std_peakloc(:,3))
scatter(none_X, std_peakloc(:,4))
errorbar(avgof_stdpeakloc(1:4), stdof_stdpeakloc(1:4), 'ob', 'linewidth', 2, 'Markersize', 10, 'MarkerFaceColor', 'b') % mean with stdev
hold off
title('Std Dev of Peak Location')
ylabel('Std Dev of Peak Location in seconds')
ax=gca;
xsize = 160
ax.YTick = [0, xsize/7, (xsize/7)*2, (xsize/7)*3, (xsize/7)*4, (xsize/7)*5, (xsize/7)*6, (xsize/7)*7];
ax.YTickLabel = {'0','1', '2', '3', '4', '5', '6', '7'};
xlabel('Stimulus Type')
ax.XTick = [1, 2, 3, 4];
ax.XTickLabel = {'Multi','Visual', 'Mech', 'None'};
fig_filename='Std Dev of Mean Peak Location RespROIs only';
saveas(gcf,fig_filename,'png');

% plot a few example ROIs (CSHL 2017 poster fig 7)
% want 1 with very low SD of multi (11 frames), average (34 frames) and
% high(59 frames)

possible_examples_M = intersect(find(std_peakloc(:,1) >= 33.5), find(std_peakloc(:,1) <=34.5))
possible_examples_L = intersect(find(std_peakloc(:,1) >= 9), find(std_peakloc(:,1) <=14))
possible_examples_H = intersect(find(std_peakloc(:,1) >= 58), find(std_peakloc(:,1) <=60))

% high var
for i = 1:length(possible_examples_H)
    figure;
    hold on 
    for j = 1:(length(goodtadpole(possible_examples_H(i)).smoothed_df_f0)*0.5)
        sweep_max = max(goodtadpole(possible_examples_H(i)).smoothed_df_f0{j})
        sweep_min = min(goodtadpole(possible_examples_H(i)).smoothed_df_f0{j})
        if (sweep_max < 0.5) && (sweep_min > -0.4)
            if goodtadpole(possible_examples_H(i)).stimorder(j) == 1
                plot(goodtadpole(possible_examples_H(i)).smoothed_df_f0{j}, 'color', [0.5 0 0.5])
            elseif goodtadpole(possible_examples_H(i)).stimorder(j) == 2
                plot(goodtadpole(possible_examples_H(i)).smoothed_df_f0{j}, 'color', 'b')
            elseif goodtadpole(possible_examples_H(i)).stimorder(j) == 3
                plot(goodtadpole(possible_examples_H(i)).smoothed_df_f0{j}, 'color', 'r')
            elseif goodtadpole(possible_examples_H(i)).stimorder(j) == 4
                plot(goodtadpole(possible_examples_H(i)).smoothed_df_f0{j}, 'color', [0.5 0.5 0.5])
            end
        end
    end
    hold off
    ax=gca;
    xsize = length(goodtadpole(possible_examples_H(i)).smoothed_df_f0{j})
    ax.XTick = [0, xsize/7, (xsize/7)*2, (xsize/7)*3, (xsize/7)*4, (xsize/7)*5, (xsize/7)*6, (xsize/7)*7];
    ax.XTickLabel = {'0','1', '2', '3', '4', '5', '6', '7'};
    title(sprintf('ROI %d all sweeps H', possible_examples_H(i)));
    xlabel('time(s)');
    ylabel('\DeltaF/F_{0}')
    fig_filename = sprintf('std_peakloc_examplesH%d', possible_examples_H(i))
    saveas(gcf,fig_filename,'epsc2');
end

%Medium var
for i = 1:length(possible_examples_M)
    figure;
    hold on 
    for j = 1:(length(goodtadpole(possible_examples_M(i)).smoothed_df_f0)*0.5)
        sweep_max = max(goodtadpole(possible_examples_M(i)).smoothed_df_f0{j})
        sweep_min = min(goodtadpole(possible_examples_M(i)).smoothed_df_f0{j})
        if length(goodtadpole(possible_examples_M(i)).smoothed_df_f0{j}) < 165
            if (sweep_max < 0.5) && (sweep_min > -0.4)
                if goodtadpole(possible_examples_M(i)).stimorder(j) == 1
                    plot(goodtadpole(possible_examples_M(i)).smoothed_df_f0{j}, 'color', [0.5 0 0.5])
                elseif goodtadpole(possible_examples_M(i)).stimorder(j) == 2
                    plot(goodtadpole(possible_examples_M(i)).smoothed_df_f0{j}, 'color', 'b')
                elseif goodtadpole(possible_examples_M(i)).stimorder(j) == 3
                    plot(goodtadpole(possible_examples_M(i)).smoothed_df_f0{j}, 'color', 'r')
                elseif goodtadpole(possible_examples_M(i)).stimorder(j) == 4
                    plot(goodtadpole(possible_examples_M(i)).smoothed_df_f0{j}, 'color', [0.5 0.5 0.5])
                end
            end
        end
    end
    hold off
    ax=gca;
    xsize = length(goodtadpole(possible_examples_M(i)).smoothed_df_f0{j})
    ax.XTick = [0, xsize/7, (xsize/7)*2, (xsize/7)*3, (xsize/7)*4, (xsize/7)*5, (xsize/7)*6, (xsize/7)*7];
    ax.XTickLabel = {'0','1', '2', '3', '4', '5', '6', '7'};
    title(sprintf('ROI %d all sweeps M', possible_examples_M(i)));
    xlabel('time(s)');
    ylabel('\DeltaF/F_{0}')
    fig_filename = sprintf('std_peakloc_examplesH%d', possible_examples_H(i))
    saveas(gcf,fig_filename,'epsc2');
end    

%low var
for i = 1:length(possible_examples_L)
    figure;
    hold on 
    for j = 1:(length(goodtadpole(possible_examples_L(i)).smoothed_df_f0)*0.5)
        sweep_max = max(goodtadpole(possible_examples_L(i)).smoothed_df_f0{j})
        sweep_min = min(goodtadpole(possible_examples_L(i)).smoothed_df_f0{j})
        if length(goodtadpole(possible_examples_L(i)).smoothed_df_f0{j}) < 165
            if (sweep_max < 0.5) && (sweep_min > -0.4)
                if goodtadpole(possible_examples_L(i)).stimorder(j) == 1
                    plot(goodtadpole(possible_examples_L(i)).smoothed_df_f0{j}, 'color', [0.5 0 0.5])
                elseif goodtadpole(possible_examples_L(i)).stimorder(j) == 2
                    plot(goodtadpole(possible_examples_L(i)).smoothed_df_f0{j}, 'color', 'b')
                elseif goodtadpole(possible_examples_L(i)).stimorder(j) == 3
                    plot(goodtadpole(possible_examples_L(i)).smoothed_df_f0{j}, 'color', 'r')
                elseif goodtadpole(possible_examples_L(i)).stimorder(j) == 4
                    plot(goodtadpole(possible_examples_L(i)).smoothed_df_f0{j}, 'color', [0.5 0.5 0.5])
                end
            end
        end
    end
    hold off
    ax=gca;
    xsize = length(goodtadpole(possible_examples_L(i)).smoothed_df_f0{j})
    ax.XTick = [0, xsize/7, (xsize/7)*2, (xsize/7)*3, (xsize/7)*4, (xsize/7)*5, (xsize/7)*6, (xsize/7)*7];
    ax.XTickLabel = {'0','1', '2', '3', '4', '5', '6', '7'};
    title(sprintf('ROI %d all sweeps L', possible_examples_L(i)));
    xlabel('time(s)');
    ylabel('\DeltaF/F_{0}')
    fig_filename = sprintf('std_peakloc_examplesL%d', possible_examples_H(i))
    saveas(gcf,fig_filename,'epsc2');
end    
    
% selected cells:
% H - 15, M - 215, L - 312
axis([5 155 -0.15 0.33])    
fig_filename = sprintf('std_peakloc_examplesM%d', 215)
saveas(gcf,fig_filename,'epsc2');    
    
% Make scatterplot of peak locations for sample ROIs
sample_rois = [15 215 312]
for k = 1:length(sample_rois)
    trials = [];
    for i = 1:4
        trials = find(goodtadpole(r).stimorder == i);
        peak_loc_ex{i} = goodtadpole(sample_rois(k)).peakloc_bytrial(1, trials);
        avg_peakloc_ex(i) = mean(peak_loc_ex{i});
        std_peakloc_ex(i) = std(peak_loc_ex{i});
    end
    multi_Xex = ones(length(peak_loc_ex{1}), 1);
    vis_Xex = ones(length(peak_loc_ex{2}), 1) + 1;
    mech_Xex = ones(length(peak_loc_ex{3}), 1) + 2;
    none_Xex = ones(length(peak_loc_ex{4}), 1) + 3;
    
    figure;
    hold on
    scatter(multi_Xex, peak_loc_ex{1}, [], [0.5 0 0.5])
    scatter(vis_Xex, peak_loc_ex{2}, [], 'b')
    scatter(mech_Xex, peak_loc_ex{3}, [], 'r')
    scatter(none_Xex, peak_loc_ex{4}, [], [0.5 0.5 0.5])
    errorbar(avg_peakloc_ex, std_peakloc_ex, 'k')
    hold off
    title(sprintf('ROI %d', sample_rois(k)))
    ylabel('Peak Location')
    ax=gca;
    ax.XTick = [1, 2, 3, 4];
    ax.XTickLabel = {'Multi', 'Vis', 'Mech', 'No stim'};
    xsize = 160
    ax.YTick = [0, xsize/7, (xsize/7)*2, (xsize/7)*3, (xsize/7)*4, (xsize/7)*5, (xsize/7)*6, (xsize/7)*7];
    ax.YTickLabel = {'0','1', '2', '3', '4', '5', '6', '7'};
    fig_filename=sprintf('std_peakloc_examples%d_scatter', sample_rois(k));
    saveas(gcf, fig_filename, 'epsc2');
    clear('peak_loc_ex')
    std_peakloc_ex_save(k,:) = std_peakloc_ex
end

% for example


        
%% Subnetworks (hypothesis 1.2, hierarchical clustering analyses)
% Tectum neurons exposed to multi-sensory conditions will give stronger
% indications of subnetworks than when exposed to a uni-sensory stimulus
% [by virtue of…]. Distinct hierarchical cluster analyses, based on each
% cell’s time-to-action potential, will be used to create subnetworks
% during multi-sensory stimulation and during uni-sensory stimulation. The
% length of each cell’s individual branch to their first cluster membership
% will be shorter in multi-sensory stimulation than in uni-sensory
% stimulation.



exps = unique([goodtadpole(:).exp]);
% multisensory peak location mean and std
for e = 1:length(exps) %cluster each exp separately
    rois = find([goodtadpole(:).exp] == exps(e))
    hierclusterMS(e).X = [std_peakloc(rois,1) avg_peakloc(rois, 1)]; 
    hierclusterMS(e).D = squareform(pdist(hierclusterMS(e).X));
    hierclusterMS(e).Z = linkage(hierclusterMS(e).D);
    hierclusterMS(e).T = cluster(hierclusterMS(e).Z, 'cutoff', 0.5);
    figure;
    dendrogram(hierclusterMS(e).Z)
    title(sprintf('Dendrogram of Multisensory Peak Location Exp %d', exps(e)))
    xlabel('ROI Number')
    ylabel('Euclidean Distance')
    fig_filename=sprintf('Dendrogram Multisensory RespROIs only exp %d', exps(e));
    saveas(gcf,fig_filename,'png');
    
end

% Visual peak location mean and std
for e = 1:length(exps) %cluster each exp separately
    rois = find([goodtadpole(:).exp] == exps(e))
    hierclusterV(e).X = [std_peakloc(rois,2) avg_peakloc(rois, 2)]; 
    hierclusterV(e).D = squareform(pdist(hierclusterV(e).X));
    hierclusterV(e).Z = linkage(hierclusterV(e).D);
    hierclusterV(e).T = cluster(hierclusterV(e).Z, 'cutoff', 0.5);
    figure;
    dendrogram(hierclusterV(e).Z)
    title(sprintf('Dendrogram of Visual Peak Location Exp %d', exps(e)))
    xlabel('ROI Number')
    ylabel('Euclidean Distance')
    fig_filename=sprintf('Dendrogram Visual RespROIs only exp %d', exps(e));
    saveas(gcf,fig_filename,'png');
    
end

% Mechanosensory peak location mean and sd
for e = 1:length(exps) %cluster each exp separately
    rois = find([goodtadpole(:).exp] == exps(e))
    hierclusterM(e).X = [std_peakloc(rois,3) avg_peakloc(rois, 3)]; 
    hierclusterM(e).D = squareform(pdist(hierclusterM(e).X));
    hierclusterM(e).Z = linkage(hierclusterM(e).D);
    hierclusterM(e).T = cluster(hierclusterM(e).Z, 'cutoff', 0.5);
    figure;
    dendrogram(hierclusterM(e).Z)
    title(sprintf('Dendrogram of Mechsensory Peak Location Exp %d', exps(e)))
    xlabel('ROI Number')
    ylabel('Euclidean Distance')
    fig_filename=sprintf('Dendrogram Mechsensory RespROIs only exp %d', exps(e));
    saveas(gcf,fig_filename,'png');
    
end

% control - no stim peak location mean and std
for e = 1:length(exps) %cluster each exp separately
    rois = find([goodtadpole(:).exp] == exps(e))
    hierclusterNS(e).X = [std_peakloc(rois,4) avg_peakloc(rois, 4)]; 
    hierclusterNS(e).D = squareform(pdist(hierclusterNS(e).X));
    hierclusterNS(e).Z = linkage(hierclusterNS(e).D);
    hierclusterNS(e).T = cluster(hierclusterNS(e).Z, 'cutoff', 0.5);
    figure;
    dendrogram(hierclusterNS(e).Z)
    title(sprintf('Dendrogram of No Stimulus Peak Location Exp %d', exps(e)))
    xlabel('ROI Number')
    ylabel('Euclidean Distance')
    fig_filename=sprintf('Dendrogram No Stimulus RespROIs only exp %d', exps(e));
    saveas(gcf,fig_filename,'png');
    
end

% just dendrograms
for e = 1:length(exps)
    figure;
    [H{e}, T{e}] = dendrogram(hierclusterNS(e).Z)
end

% Notes about Z (Linkage output)
% Z is a (m – 1)-by-3 matrix, where m is the number of observations in the
% original data. Columns 1 and 2 of Z contain cluster indices linked in
% pairs to form a binary tree. The leaf nodes are numbered from 1 to m.
% Leaf nodes are the singleton clusters from which all higher clusters are
% built. Each newly-formed cluster, corresponding to row Z(I,:), is
% assigned the index m+I. Z(I,1:2) contains the indices of the two
% component clusters that form cluster m+I. There are m-1 higher clusters
% which correspond to the interior nodes of the clustering tree. Z(I,3)
% contains the linkage distances between the two clusters merged in row
% Z(I,:). For example, suppose there are 30 initial nodes and at step 12
% cluster 5 and cluster 7 are combined. Suppose their distance at that time
% is 1.5. Then Z(12,:) will be [5, 7, 1.5]. The newly formed cluster will
% have index 12 + 30 = 42. If cluster 42 appears in a later row, it means
% the cluster created at step 12 is being combined into some larger
% cluster.

%D = pdist(X) %X = mxn data matrix, rows = observations, cols = variables
hist(hierclusterMS(5).Z(:,3))

% Statistically test to see if MS diff from US
for e = 1:length(hierclusterMS)
    [hierclusterstats(e).h_MS_V, hierclusterstats(e).p_MS_V, hierclusterstats(e).ci_MS_V hierclusterstats(e).stats_MS_V] = ttest2(hierclusterMS(e).Z(:,3), hierclusterV(e).Z(:,3))
    [hierclusterstats(e).h_MS_M, hierclusterstats(e).p_MS_M, hierclusterstats(e).ci_MS_M hierclusterstats(e).stats_MS_M] = ttest2(hierclusterMS(e).Z(:,3), hierclusterM(e).Z(:,3))
    [hierclusterstats(e).h_V_M, hierclusterstats(e).p_V_M, hierclusterstats(e).ci_V_M hierclusterstats(e).stats_V_M] = ttest2(hierclusterV(e).Z(:,3), hierclusterM(e).Z(:,3))
end
% summary: most MS-M and V-M tests have h = 1, most MS-V tests hav e h = 0.

%% Take dendrogram for e = 3 and find distances to ROIs clustered together by dendrogram algorithm
% this will only work for exps with > 30 ROIs, e = 3 has 127. 

%% Generate distance matrix first then combine based on dendrogram
for e = 1:length(exps)
    tmpidxs = find([goodtadpole(:).exp] == exps(e)); % rows of goodtadpole for that exp
    for k = 1:length(tmpidxs)
        roilocs(k, :) = goodtadpole(tmpidxs(k)).roiloc; % pull XY coordinates for ROIs in that exp
    end
    C = nchoosek((1:1:length(tmpidxs)), 2);
    for i = 1:size(C,1)
        tmpdist = sqrt( (roilocs(C(i, 1), 1) - roilocs(C(i, 2), 1))^2 + (roilocs(C(i, 1), 2) - roilocs(C(i,2), 2))^2 );
        distances{e}(i,:) = [ C(i, :) tmpdist ];
        clear('tmpdist')
    end
    clear('C', 'tmpidxs')
end

% get rois assigned to each leaf node 
for i = 1:max(T{1,3})
    Alike_ROIs{i} = find(T{1,3} == i)
end

for i = 1:length(Alike_ROIs)
    if length(Alike_ROIs{1,i}) > 1
        alike_roi_perms = nchoosek(Alike_ROIs{1,i}, 2);
        for j = 1:size(alike_roi_perms, 1)
            [possrows col] = find(distances{3}(:, 1:2) == alike_roi_perms(j,1));
            [pair_row col2] = find(distances{3}(possrows, 1:2) == alike_roi_perms(j, 2));
            Alike_ROIs{2, i}(j, :) = distances{3}(pair_row, 3);
            clear('possrows', 'col', 'col2', 'pair_row')
        end
        clear('alike_roi_perms')
    end
end

% is  the distribution of sizes for related ROIs as for the whole network?
% dump all dists into 1 vector
all_alike_dists = [];
for i = 1:size(Alike_ROIs, 2)
    if length(Alike_ROIs{2,i}) > 0
        all_alike_dists = [all_alike_dists; Alike_ROIs{2,i}];
    end
end

[h_dist, p_dist, CI_dist, stats_dist] = ttest2(distances{3}(:,3), all_alike_dists)
% p < 0.01!!!
[D PD] = allfitdist(distances{3}(:,3), 'PDF')
[D_leafs PD_leafs] = allfitdist(all_alike_dists, 'PDF')
pd = fitdist(all_alike_dists, 'Rician')
pd_all = fitdist(distances{3}(:,3), 'Rician')

figure;
hist(all_alike_dists, 20)
title('Exp 6 Resp ROIs Distances of Alike ROIs')
xlabel('Distance (um)')
ylabel('Counts')
ax=gca;
conversion = 2.26
in_um = round([0, 50, 100, 150, 200, 250] ./ conversion)
ax.XTick = [0, 50, 100, 150, 200, 250];
ax.XTickLabel = in_um;
fig_filename='hist of distances Alike ROIs';
saveas(gcf,fig_filename,'png');

figure;
hist(distances{3}(:,3), 20)
title('Exp 6 Resp ROIs Distances of All ROIs')
xlabel('Distance (um)')
ylabel('Counts')
ax=gca;
conversion = 2.26
in_um = round([0, 50, 100, 150, 200, 250] ./ conversion)
ax.XTick = [0, 50, 100, 150, 200, 250];
ax.XTickLabel = in_um;
fig_filename='hist of distances All ROIs';
saveas(gcf,fig_filename,'png');

% Map the ROI leaf nodes onto an actual tectum map
% using exp 6 only right now
%Leaf Nodes of use ( > 3 rois)
Nodes_touse = [1 2 5 7 11 21 24] 
%color_scheme = {'r', 'b', 'm', 'y', 'c', 'g', 'k'}
color_scheme = [1 0.5 0.1; 1 0 1; 0 1 0.5; 1 0.5 0; 0 1 0; 0 0.5 1; 0.5 0.5 0.5]
% roi locations (X and Y)
idxs = 172:1:298;
for i = 1:length(idxs)
    roi_locs_exp6(:,i) = [goodtadpole(idxs(i)).roiloc]
end

% primary modality
% subset of info - unimax, unimax_type,  and _multimax (peak)
primary_modality_data = [goodtadpole(182).unimax_peakavg; goodtadpole(182).unimax_stimtype; goodtadpole(182).multimax_peakavg]
for i = 1:length(primary_modality_data)
    [max_peak  max_peak_I]= max([primary_modality_data(1,i) primary_modality_data(3,i)])
    if max_peak_I == 1
        prim_mod(i) = primary_modality_data(2,i) + 1
    else
        prim_mod(i) = 1
    end
end
prim_mod_MS_all = find(prim_mod == 1)
prim_mod_M_all = find(prim_mod == 3)
prim_mod_V_all = find(prim_mod == 2)
plotted_rois = [goodtadpole(172:298).roi]
prim_mod_MS = intersect(prim_mod_MS_all, plotted_rois)
prim_mod_M = intersect(prim_mod_M_all, plotted_rois)
prim_mod_V = intersect(prim_mod_V_all, plotted_rois)
total = length(prim_mod_MS) + length(prim_mod_M) + length(prim_mod_V)
length(plotted_rois)

color_scheme_mod = [0.5 0 0.5; 0 0 1; 1 0 0]
rois_prim_mod{1} = prim_mod_MS
rois_prim_mod{2} = prim_mod_V
rois_prim_mod{3} = prim_mod_M
roi_locs_exp6(3, :) = plotted_rois
for k = 1:length(rois_prim_mod)
    for i = 1:length(rois_prim_mod{k})
        rois_prim_mod_idx{k}(i) = find(plotted_rois ==  rois_prim_mod{k}(i))
    end
end

figure;
for i = 1:length(Nodes_touse)
    hold on
    scatter(roi_locs_exp6(1, Alike_ROIs{1, Nodes_touse(i)}), roi_locs_exp6(2, Alike_ROIs{1, Nodes_touse(i)}), 50, color_scheme(i, :), 'filled')
end


for i = 1:3
    scatter(roi_locs_exp6(1, rois_prim_mod_idx{i}), roi_locs_exp6(2, rois_prim_mod_idx{i}), 50, color_scheme_mod(i,:))
end

hold off
title('ROIs by Leaf Node and Primary Modality')
fig_filename='ROIs by Leaf Node and Primary Modality Exp6';
saveas(gcf,fig_filename,'epsc2');

% Generate primary modality by node for large nodes plotted stacked bar
% graphs

for i = 1:size(Alike_ROIs, 2)
    prim_mod_bynode(i,1) = length(intersect(rois_prim_mod_idx{1}, Alike_ROIs{1, i}));
    prim_mod_bynode(i,2) = length(intersect(rois_prim_mod_idx{2}, Alike_ROIs{1, i}));
    prim_mod_bynode(i,3) = length(intersect(rois_prim_mod_idx{3}, Alike_ROIs{1, i}));
    prim_mod_bynode(i,4) = sum(prim_mod_bynode(i,1:3));
end

for i = 1:size(prim_mod_bynode, 1)
    if prim_mod_bynode(i,4) > 4
        prim_mod_bynode_prop(i, 1:3) = prim_mod_bynode(i, 1:3) ./ prim_mod_bynode(i,4);
    end
end
figure;
bar(prim_mod_bynode_prop, 'stacked')
ax=gca;
ax.XTick = [1 2 3 4 5 6];
ax.XTickLabel = {'1','2', '3', '4', '9', '12', '13'};
xlabel('Node')
ylabel('Proportion of cells')
fig_filename = 'Large Nodes Primary Modality'
saveas(gcf, fig_filename, 'epsc2')



%% Proximity analysis - for CSHL 2017 poster figure 10

% peak location/time vs roi location scatterplot for exp 6
% peak location maps onto color
% get goodtadpole(:).peakloc_avg and roilocs_exp6 into a matrix together

prox_analysis = [goodtadpole(172:298).peakloc_avg]
%prox_analysis_colormap = prox_analysis * 64/160

prox_analysis = [prox_analysis; roi_locs_exp6];
colors = [0:(1/127):1; zeros(1,128); fliplr(0:(1/127):1)]
colors_touse = colors(:, 1:127)'

colors_touse = colormap('jet')
% too short - need to map into 127 colors or turn 127 ROIs into64 colors.
max_peakloc6 = max(prox_analysis(1:3, :)')
min_peakloc6 = min(prox_analysis(1:3, :)')
% range_peakloc6 = [ 4 156 ] % range is 152

% convert the range of my peak loc values to range of 64 (so multiply by
% 64/160)
prox_analysise6MSsort_adj = floor(prox_analysise6MSsort(1,:) * (64/160))
prox_analysis_Msort_adj = floor(prox_analysis_Msort(1,:) * (64/160))
prox_analysis_Vsort_adj = floor(prox_analysis_Vsort(1,:) * (64/160))

% make figure using jet colormap. Index in by prox_analysis_Vsort_adj
figure;
scatter(prox_analysise6MSsort(2, :), prox_analysise6MSsort(3,:), 60, colors_touse(prox_analysise6MSsort_adj, :), 'filled')
colormap jet
colorbar('southoutside', 'Ticks', [0 0.5 1], 'TickLabels', {'0 sec', '3 sec', '7 sec'})
fig_filename = 'Tectum-shaped scatter by peak loc MS exp6 jet'
saveas(gcf, fig_filename, 'epsc2')


figure;
scatter(prox_analysis_Msort(2, :), prox_analysis_Msort(3,:), 60, colors_touse(prox_analysis_Msort_adj, :), 'filled')
colormap jet
colorbar('southoutside', 'Ticks', [0 0.5 1], 'TickLabels', {'0 sec', '3 sec', '7 sec'})
fig_filename = 'Tectum-shaped scatter by peak loc M exp6 jet'
saveas(gcf, fig_filename, 'epsc2')

figure;
scatter(prox_analysis_Vsort(2, :), prox_analysis_Vsort(3,:), 60, colors_touse(prox_analysis_Vsort_adj, :), 'filled')
colormap jet
colorbar('southoutside', 'Ticks', [0 0.5 1], 'TickLabels', {'0 sec', '3 sec', '7 sec'})
fig_filename = 'Tectum-shaped scatter by peak loc V exp6 jet'
saveas(gcf, fig_filename, 'epsc2')



% old figures
figure;
scatter(prox_analysise6MSsort(2, :), prox_analysise6MSsort(3,:), [], colors_touse, 'filled')
fig_filename = 'Tectum-shaped scatter by peak loc MS exp6'
saveas(gcf, fig_filename, 'epsc2')
% blue = earliest peak loc, red = latest peak loc
figure;
scatter(prox_analysis_Msort(2, :), prox_analysis_Msort(3,:), [], colors_touse, 'filled')
fig_filename = 'Tectum-shaped scatter by peak loc M exp6'
saveas(gcf, fig_filename, 'epsc2')

figure;
scatter(prox_analysis_Vsort(2, :), prox_analysis_Vsort(3,:), [], colors_touse, 'filled')
fig_filename = 'Tectum-shaped scatter by peak loc V exp6'
saveas(gcf, fig_filename, 'epsc2')

%% High/low intensity analysis - exp 22 and 23 only 
% CSHL 2017 poster figure 10
responses = sum(tadpole{6}.boolean_response_sm')
count_responses = length(find(responses))
responses23 = sum(tadpole{7}.boolean_response_sm')
count_responses = length(find(responses23))
% 10 responding ROIs total - not enough.

%% simultaneous vs offset - exp 24, tadpole 8
responses = sum(tadpole{8}.boolean_response_sm')
count_responses = length(find(responses))
% 13 responding ROIs total. not enough

%% Do "basic" analysis on just "good" rois and exps (using goodtadpole)
% Figure 6 of CSHL 2017 poster

for s = 1:4 % 4 stim types
    peaks_bytrial_all{s}.val = [];
    peaks_bytrial_all{s}.prop = [];
    tmpdata = [];
    tmpdataprop = [];
    for r = 1:length(goodtadpole) % all rois
        trials = find(goodtadpole(r).stimorder == s);
        tmpdata = [tmpdata goodtadpole(r).meanpeak_bytrial(trials)];
        tmpdataprop(r) = length(find(goodtadpole(r).meanpeak_bytrial(trials) > 0.1)) / length(trials);
        % = tmpdata2
        clear('trials')
    end
        peaks_bytrial_all{s}.val = [peaks_bytrial_all{s}.val; tmpdata];
        peaks_bytrial_all{s}.prop = tmpdataprop
        %tmpdata = [];
end

% make histograms of peak responses and proportion of trials responded too
stimtypes = {'multisensory', 'visual', 'mechanosensory', 'no stimulus'}
for s = 1:4
    type = stimtypes{s};
%     figure;
%     hist(peaks_bytrial_all{s}.val, 200)
%     title(sprintf('Peak response size all trials %d', s));
    figure;
    hist(peaks_bytrial_all{s}.prop, 20)
    title(sprintf('Proportion of trials with response > 0.1 %d', s));
    xlabel('Proportion of trials')
    ylabel('Counts')
    axis([0 1 0 80])
    clear('type')
    fig_filename = sprintf('Prop_respROIs_respond%d', s)
    saveas(gcf,fig_filename,'epsc2');
end

% are the distribution of proportions different?
% put data into a matrix for the ANOVA
for s = 1:4
    prop_resp_byroi(:,s) = peaks_bytrial_all{s}.prop;
end

[p_propresp, tbl_propresp, stats_propresp] = anova1(prop_resp_byroi)
ax=gca;
ax.XTick = [1, 2, 3, 4];
ax.XTickLabel = {'Multi', 'Vis', 'Mech', 'No stim'};
fig_filename='ANOVA of prop_responses respROIsonly';
ylabel('Proportion of trials')
saveas(gcf,fig_filename,'epsc2');

[results_propresp ,means_propresp] = multcompare(stats_propresp,'CType','bonferroni')

% Using the prop responses of good ROIs, calculate the "primary modality"
% based on highest response proportion

for r = 1:size(prop_resp_byroi,1)
    [max_val prim_mod_byprop(r)] = max(prop_resp_byroi(r, :))
end

prim_mod_byprop(2, :) = [goodtadpole(:).exp]

% Pie chart of all
for s = 1:4
    prim_mod_pie(s) = length(find(prim_mod_byprop(1,:) == s))
end

pie(prim_mod_pie, stimtypes)
fig_filename='PIE prim_mod by prop_responses respROIsonly';
saveas(gcf,fig_filename,'epsc2');

for e = 1:length(exps)
    cur = exps(e)
    cols_touse = find(prim_mod_byprop(2,:) == cur);
    for s = 1:4
        prim_mod_pie_byexp(e,s) = length(find(prim_mod_byprop(1,cols_touse) == s))
    end
    clear('cols_touse')
end

% scatter plot with mean and SD
for s = 1:4
    prim_mod_pie_byexp_percent(:,s) = prim_mod_pie_byexp(:,s) ./ sum(prim_mod_pie_byexp,2)
end
avg_prim_mod_pie_byexp_percent = mean(prim_mod_pie_byexp_percent)
std_prim_mod_pie_byexp_percent = std(prim_mod_pie_byexp_percent)
multi_X = ones(8,1)
vis_X = ones(8,1) + 1
mech_X = ones(8,1) + 2
none_X = ones(8,1) + 3

figure;
hold on
% all percentages as circles
scatter(multi_X, prim_mod_pie_byexp_percent(:,1))
scatter(vis_X, prim_mod_pie_byexp_percent(:,2))
scatter(mech_X, prim_mod_pie_byexp_percent(:,3))
scatter(none_X, prim_mod_pie_byexp_percent(:,4))
errorbar(avg_prim_mod_pie_byexp_percent, std_prim_mod_pie_byexp_percent, 'k')
hold off
ylabel('Proportion of cells')
ax=gca;
ax.XTick = [1, 2, 3, 4];
ax.XTickLabel = {'Multi', 'Vis', 'Mech', 'No stim'};
fig_filename='Percent prim_mod_byexp';

saveas(gcf,fig_filename,'epsc2');

% Change primary modality pie chart + scatterplot into stacked bar graph by
% exp
bar(prim_mod_pie_byexp_percent, 'stacked')
ylabel('Proportion of cells')
fig_filename = 'Primary modality all exps stacked bar'
saveas(gcf, fig_filename, 'epsc2')


% for example

figure;

title('Std Dev of Peak Location')
ylabel('Std Dev of Peak Location in seconds')
ax=gca;
xsize = 160
ax.YTick = [0, xsize/7, (xsize/7)*2, (xsize/7)*3, (xsize/7)*4, (xsize/7)*5, (xsize/7)*6, (xsize/7)*7];
ax.YTickLabel = {'0','1', '2', '3', '4', '5', '6', '7'};
xlabel('Stimulus Type')
ax.XTick = [1, 2, 3, 4];
ax.XTickLabel = {'Multi','Visual', 'Mech', 'None'};
fig_filename='Std Dev of Mean Peak Location RespROIs only';
saveas(gcf,fig_filename,'png');



%% Map dendrogram onto physical distance (DOESN'T WORK AND IDK WHY)
% Take the first connection for each ROI and get the physical Euclidean
% distance between the ROIs. If multiple ROIs, take average of distances.
% Determine if physical distance correlates with dendrogram distance. Is
% this affected by significance of  ttest above?

for e = 1:length(exps)
    tmpidxs = find([goodtadpole(:).exp] == exps(e)); % rows of goodtadpole for that exp
    for k = 1:length(tmpidxs)
        roilocs(k, :) = goodtadpole(tmpidxs(k)).roiloc; % pull XY coordinates for ROIs in that exp
    end
    len = size(hierclusterMS(e).X, 1) % number of rois
    for i = 1:len
        roi1 = hierclusterMS(e).Z(i, 1)
        roi2 = hierclusterMS(e).Z(i, 2)
        if (roi1 < len) && (roi2 < len)
            hierclusterMS(e).Phys_dist(i) = sqrt( (roilocs(roi1,1) - roilocs(roi2, 1))^2 + (roilocs(roi1,2) - roilocs(roi2, 2))^2 );
        elseif (roi1 > len) && (roi2 < len)
            roi1a = hierclusterMS(e).Z((roi1-len), 1:2);
            if roi1a < len % if both are less than len
                dist1 = sqrt( (roilocs(roi1a(1),1) - roilocs(roi2, 1))^2 + (roilocs(roi1a(1),2) - roilocs(roi2, 2))^2 );
                dist2 = sqrt( (roilocs(roi1a(2),1) - roilocs(roi2, 1))^2 + (roilocs(roi1a(2),2) - roilocs(roi2, 2))^2 );
                hierclusterMS(e).Phys_dist(i) = mean([dist1 dist2]);
            else
                largeidx = find(roi1a > len) % idx of too big roi num
                smallidx = find(roi1a < len) % idx of normal roi num
                roi1b = hierclusterMS(e).Z((roi1a(largeidx)-len), 1:2);
                dist1 = sqrt( (roilocs(roi1b(1),1) - roilocs(roi2, 1))^2 + (roilocs(roi1b(1),2) - roilocs(roi2, 2))^2 );
                dist2 = sqrt( (roilocs(roi1b(2),1) - roilocs(roi2, 1))^2 + (roilocs(roi1b(2),2) - roilocs(roi2, 2))^2 );
                avg_roi1b = mean([dist1 dist2]); % average the 2 distances that make up roi1a.
                dist3 = sqrt( (roilocs(roi1a(smallidx),1) - roilocs(roi2, 1))^2 + (roilocs(roi1a(smallidx),2) - roilocs(roi2, 2))^2 );
                hierclusterMS(e).Phys_dist(i) = mean([avg_roi1b dist3]); % average the distances that make up roi1
            end
        elseif (roi1 < len) && (roi2 > len)
            roi2a = hierclusterMS(e).Z((roi2-len), 1:2);
            if roi2a < len % if both are less than len
                dist1 = sqrt( (roilocs(roi1,1) - roilocs(roi2a(1), 1))^2 + (roilocs(roi1,2) - roilocs(roi2a(1), 2))^2 );
                dist2 = sqrt( (roilocs(roi1,1) - roilocs(roi2a(2), 1))^2 + (roilocs(roi1,2) - roilocs(roi2a(2), 2))^2 );
                hierclusterMS(e).Phys_dist(i) = mean([dist1 dist2]);
             else
                largeidx = find(roi2a > len) % idx of too big roi num
                smallidx = find(roi2a < len) % idx of normal roi num
                roi2b = hierclusterMS(e).Z((roi2a(largeidx)-len), 1:2);
                dist1 = sqrt( (roilocs(roi1,1) - roilocs(roi2b(1), 1))^2 + (roilocs(roi1,2) - roilocs(roi2b(1), 2))^2 );
                dist2 = sqrt( (roilocs(roi1,1) - roilocs(roi2b(2), 1))^2 + (roilocs(roi1,2) - roilocs(roi2b(2), 2))^2 );
                avg_roi2b = mean([dist1 dist2]); % average the 2 distances that make up roi1a.
                dist3 = sqrt( (roilocs(roi1,1) - roilocs(roi2a(smallidx), 1))^2 + (roilocs(roi1,2) - roilocs(roi2a(smallidx), 2))^2 );
                hierclusterMS(e).Phys_dist(i) = mean([avg_roi2b dist3]); % average the distances that make up roi1    
            end
       elseif (roi1 > len) && (roi2 > len)
            roi1a = hierclusterMS(e).Z((roi1-len), 1:2);
            roi2a = hierclusterMS(e).Z((roi2-len), 1:2);
            if (roi1a < len) && (roi2a < len) % both rois nested in roi1 and roi2 are real rois
                dist1 = sqrt( (roilocs(roi1a(1),1) - roilocs(roi2a(1), 1))^2 + (roilocs(roi1a(1),2) - roilocs(roi2a(1), 2))^2 );
                dist2 = sqrt( (roilocs(roi1a(2),1) - roilocs(roi2a(1), 1))^2 + (roilocs(roi1a(2),2) - roilocs(roi2a(1), 2))^2 );
                dist3 = sqrt( (roilocs(roi1a(1),1) - roilocs(roi2a(2), 1))^2 + (roilocs(roi1a(1),2) - roilocs(roi2a(2), 2))^2 );
                dist4 = sqrt( (roilocs(roi1a(2),1) - roilocs(roi2a(2), 1))^2 + (roilocs(roi1a(2),2) - roilocs(roi2a(2), 2))^2 );
                hierclusterMS(e).Phys_dist(i) = mean([dist1 dist2 dist3 dist4]);
            
            elseif roi1a > len % roi2a are real rois
                largeidx1 = find(roi1a > len) % idx of too big roi1 num
                smallidx1 = find(roi1a < len) % idx of normal roi1 num
                % 1b1-2a1, 1b1-2a2, 1b2-2a1, 1b2-2a2
                roi1b = hierclusterMS(e).Z((roi1a(largeidx1)-len), 1:2);
                % roi1b 1, 2 against 2a 1,2
                dist1 = sqrt( (roilocs(roi1b(1),1) - roilocs(roi2a(1), 1))^2 + (roilocs(roi1b(1),2) - roilocs(roi2a(1), 2))^2 );
                dist2 = sqrt( (roilocs(roi1b(2),1) - roilocs(roi2a(1), 1))^2 + (roilocs(roi1b(2),2) - roilocs(roi2a(1), 2))^2 );
                dist3 = sqrt( (roilocs(roi1b(1),1) - roilocs(roi2a(2), 1))^2 + (roilocs(roi1b(1),2) - roilocs(roi2a(2), 2))^2 );
                dist4 = sqrt( (roilocs(roi1b(2),1) - roilocs(roi2a(2), 1))^2 + (roilocs(roi1b(2),2) - roilocs(roi2a(2), 2))^2 );
                avg_roi1b_roi2a = mean([dist1 dist2 dist3 dist4]);
                % roi1a(smallidx) against 2a 1, 2
                dist5 = sqrt( (roilocs(roi1a(smallidx1),1) - roilocs(roi2a(1), 1))^2 + (roilocs(roi1a(smallidx1),2) - roilocs(roi2a(1), 2))^2 );
                dist6 = sqrt( (roilocs(roi1a(smallidx1),1) - roilocs(roi2a(2), 1))^2 + (roilocs(roi1a(smallidx1),2) - roilocs(roi2a(2), 2))^2 );
                avg_roi1a_roi2b = mean([dist5 dist6]);
                
            elseif roi2a > len % roi1a are real rois   
                largeidx2 = find(roi2a > len) % idx of too big roi2 num
                smallidx2 = find(roi2a < len) % idx of normal roi2 num                
            end
        end
    end
    clear('tmpidxs', 'roilocs')
end
end



% Find distance for dendrogram
for e = 1:length(exps)
    len = size(hierclusterMS(e).X, 1)
    for i = 1:len % over all rois. i is roi1, then look for other rois. 
        [row col roi1] = find(hierclusterMS(e).Z(:, 1:2) == i)
        if col == 1 
            roi2 = hierclusterMS(e).Z(row, 2)
        elseif col == 2
            roi2 = hierclusterMS(e).Z(row, 1)
        end
        smaller_roi = min([i, roi2])
        larger_roi = max([i, roi2])
        if larger_roi > len
            larger_roiA = hierclusterMS(e).Z((larger_roi-len), 1:2);
            %if larger_roiA < len
                possrows1 = find(distances{e}(:, 1) == smaller_roi)
                distrow1 = find(distances{e}(possrows1, 2) == larger_roiA(1))
                possrows2 = find(distances{e}(:, 1) == smaller_roi)
                distrow2 = find(distances{e}(possrows2, 2) == larger_roiA(2))
                firstdistance{e}(i, :) = [ smaller_roi, larger_roi, mean([distances{e}(distrow1, 3), distances{e}(distrow2, 3)]) ]; 
           %else 
            %    larger_roiAB
        else    
            possrows = find(distances{e}(:, 1) == smaller_roi)
            distrow = find(distances{e}(possrows, 2) == larger_roi)
            firstdistance{e}(i, :) = [ smaller_roi, larger_roi, distances{e}(distrow, 3) ]; 
        end
    end
end

        
        
        
        
        
        
        
        
        
        roi1 = hierclusterMS(e).Z(i, 1)
        roi2 = hierclusterMS(e).Z(i, 2)
        if (roi1 < len) && (roi2 < len)
            [row1 col1] = find(distances{e}(:, 1:2) == roi1)
            [row2 col2] = find(distances{e}(:, 1:2) == roi2)
            rowid = intersect(row1, row2)
            hiercluster(e).distances(i) = distances{e}(rowid, 3);
        elseif (roi1 > len)  
            roi1a = hierclusterMS(e).Z((roi1-len), 1:2);
            if roi1a < len
                if roi2 < len
                    [row1a1 col1a1] = find(distances{e}(:, 1:2) == roi1a(1))
                    [row1a2 col1a2] = find(distances{e}(:, 1:2) == roi1a(2))
                    [row2 col2] = find(distances{e}(:, 1:2) == roi2)
                    rowid1a1 = intersect(row1a1, row2)
                    rowid1a2 = intersect(row1a2, row2)
                    hiercluster(e).distances(i) = mean([distances{e}(rowid1a1, 3), distances{e}(rowid1a2, 3)]);
                elseif roi2 > len
                    roi2a = hierclusterMS(e).Z((roi2-len), 1:2);
                    [row1a1 col1a1] = find(distances{e}(:, 1:2) == roi1a(1))
                    [row1a2 col1a2] = find(distances{e}(:, 1:2) == roi1a(2))
                    [row2a1 col2a1] = find(distances{e}(:, 1:2) == roi2a(1))
                    [row2a2 col2a2] = find(distances{e}(:, 1:2) == roi2a(2))
                    rowid1a1_2a1 = intersect(row1a1, row2a1)
                    rowid1a2_2a1 = intersect(row1a2, row2a1)
                    rowid1a1_2a2 = intersect(row1a1, row2a2)
                    rowid1a2_2a2 = intersect(row1a2, row2a2)
                    hiercluster(e).distances(i) = mean([distances{e}(rowid1a1_2a1, 3), distances{e}(rowid1a2_2a1, 3), distances{e}(rowid1a1_2a2, 3), distances{e}(rowid1a2_2a2, 3)]);
                end
            elseif roi1a > len
                largeidx = find(roi1a > len) % idx of too big roi1 num
                smallidx = find(roi1a < len) % idx of normal roi1 num
                roi1b = hierclusterMS(e).Z((roi1a(largeidx-len), 1:2);
                if roi2 < len
                    [row1a1 col1a1] = find(distances{e}(:, 1:2) == roi1a(smallidx))
                    [row1a2b1 col1a2b1] = find(distances{e}(:, 1:2) == roi1b(1))
                    [row1a2b2 col1a2b2] = find(distances{e}(:, 1:2) == roi1b(2))
                    [row2 col2] = find(distances{e}(:, 1:2) == roi2)
                    rowid1a1 = intersect(row1a1, row2)
                    rowid1a2 = intersect(row1a2, row2)
                    hiercluster(e).distances(i) = mean([distances{e}(rowid1a1, 3), distances{e}(rowid1a2, 3)]);
                end    
            end
        end

        
        
        for i = 1:length(hierclusterMS)
            roinum(i) = size(hierclusterMS(i).X, 1)
            largestclusterval(i) = max(max(hierclusterMS(i).Z(:, 1:2)))
        end
        largestclusterval./roinum

                
