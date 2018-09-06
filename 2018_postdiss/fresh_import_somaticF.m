%% Get all raw data into 1 struct (not a cell array of structs as I have made before)

% Collected raw data for tads 2-11 from F:\Calcium_Imaging_Analysis\analyzed_compiled\Smoothed_analysis
%% Import tadpole struct for all experiments
myFolder = 'D:\\Torrey_calcium_imaging\2018August\raw_data' %'F:/Calcium_Imaging_Analysis/analyzed_compiled/Smoothed_analysis/'; % May need to correct this.
%mkdir([myFolder 'figures']);
if ~isdir(myFolder)
	errorMessage = sprintf('Error: The following folder does not exist:\n%s', myFolder);
	uiwait(warndlg(errorMessage));
	return;
end
filePattern = fullfile(myFolder, 'tad*.mat');
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
    %%%%% old
	%tadpole{1,k} = matData.tadpole % If you get to here, tadpole existed in the file.
    %tadpole{1,k}.matFilename = matFilename
    
    %%%%%% new
    allData(k).expnum = matData.tadpole.expnum;
    %allData(k).expdate = matData.tadpole.expdate;
    allData(k).trial_length = matData.tadpole.trial_length;
    allData(k).numtrialblocks = matData.tadpole.numtrialblocks;
    allData(k).num_trials = matData.tadpole.num_trials;
    %allData(k).blockIDs = matData.tadpole.blockids;
    allData(k).somaticF = matData.tadpole.somaticF;
    allData(k).somaticROICenters = matData.tadpole.somaticROICenters;
    allData(k).neuropilF = matData.tadpole.neuropilF;
    allData(k).stimorder = matData.tadpole.stimorder;
    allData(k).trial_splitS = matData.tadpole.trial_splitS;
    allData(k).trial_splitN = matData.tadpole.trial_splitN;
    allData(k).background = matData.tadpole.background;
    allData(k).backgroundsubS = matData.tadpole.backgroundsubS;
    allData(k).signal = matData.tadpole.signal;
    allData(k).df_f0 = matData.tadpole.df_f0;
    %allData(k).smoothed = matData.tadpole.smoothed;
    
end

%% %%%%%%%%% WHAT FIELDS DO I HAVE?
for k = 1:length(matFiles)
	matFilename = fullfile(myFolder, matFiles(k).name)
	matData = load(matFilename, 'tadpole'); % Retrieves a structure.

    if isfield(matData.tadpole, 'expnum')
        what_fields(k,1) = 1;
    end
    if isfield(matData.tadpole, 'expdate')
        what_fields(k,2) = 1;
    end
    if isfield(matData.tadpole, 'trial_length')
        what_fields(k,3) = 1;
    end    
    if isfield(matData.tadpole, 'numtrialblocks')
        what_fields(k,4) = 1;
    end
        if isfield(matData.tadpole, 'num_trials')
        what_fields(k,5) = 1;
        end
        if isfield(matData.tadpole, 'blockids')
        what_fields(k,6) = 1;
        end
        if isfield(matData.tadpole, 'somaticF')
        what_fields(k,7) = 1;
        end
        if isfield(matData.tadpole, 'somaticROICenters')
        what_fields(k,8) = 1;
        end
        if isfield(matData.tadpole, 'neuropilF')
        what_fields(k,9) = 1;
        end
        if isfield(matData.tadpole, 'stimorder')
        what_fields(k,10) = 1;
        end
        if isfield(matData.tadpole, 'trial_splitS')
        what_fields(k,11) = 1;
        end
        if isfield(matData.tadpole, 'trial_splitN')
        what_fields(k,12) = 1;
        end    
        if isfield(matData.tadpole, 'background')
        what_fields(k,13) = 1;
        end
        if isfield(matData.tadpole, 'backgroundsubS')
        what_fields(k,14) = 1;
        end
        if isfield(matData.tadpole, 'signal')
        what_fields(k,15) = 1;
        end
        if isfield(matData.tadpole, 'df_f0')
        what_fields(k,16) = 1;
        end
        if isfield(matData.tadpole, 'smoothed')
        what_fields(k,17) = 1;
        end
%     %%%%%% new
%     allData(k).expnum = matData.tadpole.expnum;
%     %allData(k).expdate = matData.tadpole.expdate;
%     allData(k).trial_length = matData.tadpole.trial_length;
%     allData(k).numtrialblocks = matData.tadpole.numtrialblocks;
%     allData(k).num_trials = matData.tadpole.num_trials;
%     %allData(k).blockIDs = matData.tadpole.blockids;
%     allData(k).somaticF = matData.tadpole.somaticF;
%     allData(k).somaticROICenters = matData.tadpole.somaticROICenters;
%     allData(k).neuropilF = matData.tadpole.neuropilF;
%     allData(k).stimorder = matData.tadpole.stimorder;
%     allData(k).trial_splitS = matData.tadpole.trial_splitS;
%     allData(k).trial_splitN = matData.tadpole.trial_splitN;
%     allData(k).background = matData.tadpole.background;
%     allData(k).backgroundsubS = matData.tadpole.backgroundsubS;
%     allData(k).signal = matData.tadpole.signal;
%     allData(k).df_f0 = matData.tadpole.df_f0;
%     %allData(k).smoothed = matData.tadpole.smoothed;
    
end
what_fields(:,18) = [allData(:).expnum];
what_fields(31,:) = sum(what_fields(1:30,:));

%% Import everything that exists for a given tadpole:

for k = 1:length(matFiles)
	matFilename = fullfile(myFolder, matFiles(k).name)
	matData = load(matFilename, 'tadpole'); % Retrieves a structure.

    if isfield(matData.tadpole, 'expnum')
        allData(k).expnum = matData.tadpole.expnum;
    end
    if isfield(matData.tadpole, 'expdate')
        allData(k).expdate = matData.tadpole.expdate;
    end
    if isfield(matData.tadpole, 'trial_length')
        allData(k).trial_length = matData.tadpole.trial_length;
    end    
    if isfield(matData.tadpole, 'numtrialblocks')
        allData(k).numtrialblocks = matData.tadpole.numtrialblocks;
    end
    if isfield(matData.tadpole, 'num_trials')
        allData(k).num_trials = matData.tadpole.num_trials;
    end
    if isfield(matData.tadpole, 'blockids')
        allData(k).blockIDs = matData.tadpole.blockids;
    end
    if isfield(matData.tadpole, 'somaticF')
        allData(k).somaticF = matData.tadpole.somaticF;
    end
    if isfield(matData.tadpole, 'somaticROICenters')
        allData(k).somaticROICenters = matData.tadpole.somaticROICenters;
    end
    if isfield(matData.tadpole, 'neuropilF')
        allData(k).neuropilF = matData.tadpole.neuropilF;
    end
    if isfield(matData.tadpole, 'stimorder')
        allData(k).stimorder = matData.tadpole.stimorder;
    end
    if isfield(matData.tadpole, 'trial_splitS')
        allData(k).trial_splitS = matData.tadpole.trial_splitS;
    end
    if isfield(matData.tadpole, 'trial_splitN')
        allData(k).trial_splitN = matData.tadpole.trial_splitN;
    end    
    if isfield(matData.tadpole, 'background')
        allData(k).background = matData.tadpole.background;
    end
    if isfield(matData.tadpole, 'backgroundsubS')
        allData(k).backgroundsubS = matData.tadpole.backgroundsubS;
    end
    if isfield(matData.tadpole, 'signal')
        allData(k).signal = matData.tadpole.signal;
    end
    if isfield(matData.tadpole, 'df_f0')
        allData(k).df_f0 = matData.tadpole.df_f0;
    end
    if isfield(matData.tadpole, 'smoothed')
        allData(k).smoothed = matData.tadpole.smoothed;
    end
end

%% Now that I have all (or at least minimum to reconstruct everything), fill in any missing data points
for k = 1:length(allData)
    if isempty(allData(k).smoothed)
        for i = 1:size(allData(k).df_f0, 1)
            for j = 1:size(allData(k).df_f0, 2)
                allData(k).smoothed{i,j} = smooth(allData(k).df_f0{i,j}(:,:), 8, 'moving');
            end
        end
    end
end

%% Now that I have raw mean fluorescence through smoothed, make plots 

% plot all ROIs for each trial, with each type in a separate subplot
k = 10
tr = 10
for k = 1:length(allData) % each tadpole)
   roi_ct = size(allData(k).trial_splitS, 1)
   for tr = 1:size(allData(k).trial_splitS, 2) %each trial in tadpole
        figure;
        % top - raw mean fluorescence (trial_splitS)
        subplot(4,1,1)
        hold on
        for j = 1:roi_ct
            plot(allData(k).trial_splitS{j, tr})
        end
        hold off
        set(gca,'xtick',[])
        ylabel('raw')
        title(sprintf('All Signal Steps tad %d (k=%d) trial %d', allData(k).expnum, k, tr));
        
        % second - background and neuropil subtracted (signal)
        subplot(4,1,2)
        hold on
        for j = 1:roi_ct
            plot(allData(k).signal{j, tr})
        end
        hold off
        set(gca,'xtick',[])
        ylabel('signal')
        
        % third - df/f0
        subplot(4,1,3)
        hold on
        for j = 1:roi_ct
            plot(allData(k).df_f0{j, tr})
        end
        hold off
        set(gca,'xtick',[])
        ylabel('df/f_{0}')
        
        % bottom - smoothed df/f0 (8 point moving average)
        subplot(4,1,4)
        hold on
        for j = 1:roi_ct
            plot(allData(k).smoothed{j, tr})
        end
        hold off
        ylabel('sm df/f_{0}')
        ax=gca;
        xsize = length(allData(k).smoothed{j, tr});
        ax.XTick = [0, xsize/7, (xsize/7)*2, (xsize/7)*3, (xsize/7)*4, (xsize/7)*5, (xsize/7)*6, (xsize/7)*7];
        ax.XTickLabel = {'0','1', '2', '3', '4', '5', '6', '7'};
        xlabel('time (s)');
        
        fig_filename=sprintf('all_signal_steps_tad%d_k%d_trial%d', allData(k).expnum, k, tr);
        saveas(gcf,fig_filename,'png');
        saveas(gcf, fig_filename, 'epsc2');
        close;
        clear('fig_filename')
   end
end

%% Ok, so the spurious values are appearing at the calculation of df/f and being passed through to smoothed df/f

% Therefore - need to figure out the source. Why are some traces going
% screwy?

% first thought - f0 values very close to 0 would cause spurious large values. 
% dummy test
test_df = [1 2 3 4 5 6 7 8 9 10]
test_f0 = [0.001 0.01 0.1 1 10]

for d = 1:length(test_df)
    for f = 1:length(test_f0)
        test_df_f0 (d,f) = (test_df(d) - test_f0(f))/test_f0(f);
    end
end

% result - yes. f0 = 0.001 gives vals in thousands. 

% Could this be just the background ROI being included?
% test tad 6 trial 4 (k = 27)
k = 27
tr = 4
roi_ct = size(allData(k).trial_splitS, 1)-1
figure;
        % top - raw mean fluorescence (trial_splitS)
        subplot(4,1,1)
        hold on
        for j = 1:roi_ct
            plot(allData(k).trial_splitS{j, tr})
        end
        hold off
        set(gca,'xtick',[])
        ylabel('raw')
        title(sprintf('All Signal Steps tad %d (k=%d) trial %d', allData(k).expnum, k, tr));
        
        % second - background and neuropil subtracted (signal)
        subplot(4,1,2)
        hold on
        for j = 1:roi_ct
            plot(allData(k).signal{j, tr})
        end
        hold off
        set(gca,'xtick',[])
        ylabel('signal')
        
        % third - df/f0
        subplot(4,1,3)
        hold on
        for j = 1:roi_ct
            plot(allData(k).df_f0{j, tr})
        end
        hold off
        set(gca,'xtick',[])
        ylabel('df/f_{0}')
        
        % bottom - smoothed df/f0 (8 point moving average)
        subplot(4,1,4)
        hold on
        for j = 1:roi_ct
            plot(allData(k).smoothed{j, tr})
        end
        hold off
        ylabel('sm df/f_{0}')
        ax=gca;
        xsize = length(allData(k).smoothed{j, tr});
        ax.XTick = [0, xsize/7, (xsize/7)*2, (xsize/7)*3, (xsize/7)*4, (xsize/7)*5, (xsize/7)*6, (xsize/7)*7];
        ax.XTickLabel = {'0','1', '2', '3', '4', '5', '6', '7'};
        xlabel('time (s)');
% NOPE (or, at least not entirely)

% so what ROI is creating the problem here?
for r = 1:size(allData(k).df_f0,1)
    tmp = find(allData(k).df_f0{r,tr} > 1);
    if sum(tmp) > 0
        spurious_trace(r) = 1;
    end
end
rois = find(spurious_trace)
% roi are 3 random ROIs 

% plot those spurious traces by themselves
rois = [51 56 108];

figure;
% top - raw mean fluorescence (trial_splitS)
subplot(4,1,1)
hold on
for j = 1:length(rois)
    plot(allData(k).trial_splitS{rois(j), tr})
end
hold off
set(gca,'xtick',[])
ylabel('raw')
title(sprintf('All Signal Steps tad %d (k=%d) trial %d', allData(k).expnum, k, tr));

% second - background and neuropil subtracted (signal)
subplot(4,1,2)
hold on
for j = 1:length(rois)
    plot(allData(k).signal{rois(j), tr})
end
hold off
set(gca,'xtick',[])
ylabel('signal')

% third - df/f0
subplot(4,1,3)
hold on
for j = 1:length(rois)
    plot(allData(k).df_f0{rois(j), tr})
end
hold off
set(gca,'xtick',[])
ylabel('df/f_{0}')

% bottom - smoothed df/f0 (8 point moving average)
subplot(4,1,4)
hold on
for j = 1:length(rois)
    plot(allData(k).smoothed{rois(j), tr})
end
hold off
ylabel('sm df/f_{0}')
ax=gca;
xsize = length(allData(k).smoothed{j, tr});
ax.XTick = [0, xsize/7, (xsize/7)*2, (xsize/7)*3, (xsize/7)*4, (xsize/7)*5, (xsize/7)*6, (xsize/7)*7];
ax.XTickLabel = {'0','1', '2', '3', '4', '5', '6', '7'};
xlabel('time (s)');

% Yes, ROIs 51, 56 and 108 have "signal" values very close to 0. 

% calculate all F0 values 
for k = 1:length(allData) % all tads
    for r = 1:size(allData(k).signal,1) %all rois
        for tr = 1:size(allData(k).signal,2) % all trials
            allData(k).f0(r, tr) = mean(allData(k).signal{r,tr}(1:10));
        end
    end
end




% % Are traces with f0~=0 likely to have very large or very small df peak
% % values?
% 
% 
% for i = 1:size(signal,1)
%     for j = 1:size(signal,2)
%         F0 = mean(signal{i,j}(1:stimulus_frame));
%         df_f0{i,j} = (signal{i,j} - F0)/F0;
%     end
% end

%% Replot all data by ROI - same as above but put all trials for an ROI on 1 plot, rather than all ROIs for a given trial
%k = 10
%r = 10
for k = 1:length(allData) % each tadpole)
   trial_ct = size(allData(k).trial_splitS, 2)
   for r = 1:size(allData(k).trial_splitS, 1) %each trial in tadpole
        figure;
        % top - raw mean fluorescence (trial_splitS)
        subplot(4,1,1)
        hold on
        for j = 1:trial_ct
            plot(allData(k).trial_splitS{r, j})
            min_val(j) = min(allData(k).trial_splitS{r, j});
            max_val(j) = max(allData(k).trial_splitS{r, j});
        end
        hold off
        set(gca,'xtick',[])
        ylabel('raw')
        ylim([min(min_val) max(max_val)])
        clear('min_val', 'max_val')
        title(sprintf('All Signal Steps tad %d (k=%d) roi %d', allData(k).expnum, k, r));
        
        % second - background and neuropil subtracted (signal)
        subplot(4,1,2)
        hold on
        for j = 1:trial_ct
            plot(allData(k).signal{r, j})
            min_val(j) = min(allData(k).signal{r, j});
            max_val(j) = max(allData(k).signal{r, j});
        end
        hold off
        set(gca,'xtick',[])
        ylabel('signal')
        ylim([min(min_val) max(max_val)])
        clear('min_val', 'max_val')
        
        % third - df/f0
        subplot(4,1,3)
        hold on
        for j = 1:trial_ct
            plot(allData(k).df_f0{r, j})
            min_val(j) = min(allData(k).df_f0{r, j});
            max_val(j) = max(allData(k).df_f0{r, j});
        end
        hold off
        set(gca,'xtick',[])
        ylabel('df/f_{0}')
        ylim([min(min_val) max(max_val)])
                
        % bottom - smoothed df/f0 (8 point moving average)
        subplot(4,1,4)
        hold on
        for j = 1:trial_ct
            plot(allData(k).smoothed{r, j})
        end
        hold off
        ylabel('sm df/f_{0}')
        ylim([min(min_val) max(max_val)]) %same as original df/f
        clear('min_val', 'max_val')
        
        % display seconds on x axis rather than frames
        ax=gca;
        xsize = length(allData(k).smoothed{r, j});
        ax.XTick = [0, xsize/7, (xsize/7)*2, (xsize/7)*3, (xsize/7)*4, (xsize/7)*5, (xsize/7)*6, (xsize/7)*7];
        ax.XTickLabel = {'0','1', '2', '3', '4', '5', '6', '7'};
        xlabel('time (s)');
        
        % save 
        fig_filename=sprintf('all_signal_steps_tad%d_k%d_roi%d', allData(k).expnum, k, r);
        saveas(gcf,fig_filename,'png');
        saveas(gcf, fig_filename, 'epsc2');
        close;
        clear('fig_filename')
   end
end

%% Check smoothing algorithm 
% Chris Moore said it's technically more correct to run the smoothing in
% both directions. I only ran it forwards in allData(k).smoothed

for k = 1:length(allData)
    for i = 1:size(allData(k).smoothed, 1)
        for j = 1:size(allData(k).smoothed, 2)
            allData(k).smoothed_bd{i, j} = flipud(smooth(flipud(allData(k).smoothed{i,j}(:,:)), 8, 'moving'));
        end
    end
end

% %test plot
% figure;
% hold on
% %for i = 1:size(allData(k).df_f0,1)
% plot(allData(k).df_f0{i, j}, 'k')
% plot(allData(k).smoothed{i, j}, 'g')
% plot(allData(k).smoothed_bd{i, j}, 'r')
% hold off

%% Plot df/f, smoothed and smoothed_bd (both directions)

for k = 1:length(allData) % each tadpole)
   roi_ct = size(allData(k).df_f0, 1)
   for tr = 1:size(allData(k).df_f0, 2) %each trial in tadpole
        figure;
        % top - df_f0 
        subplot(3,1,1)
        hold on
        for j = 1:roi_ct
            plot(allData(k).df_f0{j, tr})
            min_val(j) = min(allData(k).df_f0{j, tr});
            max_val(j) = max(allData(k).df_f0{j, tr});
        end
        hold off
        set(gca,'xtick',[])
        ylim([min(min_val) max(max_val)])
        %clear('min_val', 'max_val')
        ylabel('df/f_{0}')
        title(sprintf('Smoothed and sm_bd tad %d (k=%d) trial %d', allData(k).expnum, k, tr));
        
        % second - background and neuropil subtracted (signal)
        subplot(3,1,2)
        hold on
        for j = 1:roi_ct
            plot(allData(k).smoothed{j, tr});
        end
        hold off
        set(gca,'xtick',[])
        ylabel('sm df/f_{0}')
        ylim([min(min_val) max(max_val)])
        
        % third - df/f0
        subplot(3,1,3)
        hold on
        for j = 1:roi_ct
            plot(allData(k).smoothed_bd{j, tr})
        end
        hold off
        %set(gca,'xtick',[])
        ylabel('bd sm df/f_{0}')
        ylim([min(min_val) max(max_val)])
        clear('min_val', 'max_val')
        
        ax=gca;
        xsize = length(allData(k).smoothed{j, tr});
        ax.XTick = [0, xsize/7, (xsize/7)*2, (xsize/7)*3, (xsize/7)*4, (xsize/7)*5, (xsize/7)*6, (xsize/7)*7];
        ax.XTickLabel = {'0','1', '2', '3', '4', '5', '6', '7'};
        xlabel('time (s)');
        
        fig_filename=sprintf('each_smoothed_tad%d_k%d_trial%d', allData(k).expnum, k, tr);
        saveas(gcf,fig_filename,'png');
        saveas(gcf, fig_filename, 'epsc2');
        close;
        clear('fig_filename')
   end
end

%%%% Looking at the resulting plots, smoothing both directions just dilutes
%%%% the signal - does not take care of the weird looking first few and
%%%% last few points. 

%% Ok, so let's look at the raw df/f (no smoothing) instead
%%% Going to require better quality traces to recover useful info.

%% First, should I bleach correct to reduce the effect of fluorophore decay over a trial

% Assess spontaneous (4) trials to determine rate of bleaching for each ROI
% in each tad

% going to assume linear decay rate over this short time frame

for k = 1:length(allData) % over all tads
    spont_tr = find(allData(k).stimorder == 4);
    for t = 1:length(spont_tr) %trials
        len = 1:1:length(allData(k).signal{1,t});
        for r = 1:size(allData(k).signal,1) % rois
            allData(k).decay(r,t,:) = polyfit(allData(k).signal{r,t}, len, 1);
        end
    end
end

% Is there a change in the decay constant over trials?
for k = 1:length(allData)
    figure;
    hold on
    plot(allData(k).decay(:,:,1)')
    plot([1:1:size(allData(k).decay,2)], mean(allData(k).decay(:,:,1), 1), 'k', 'LineWidth', 2)
    xlim([1, size(allData(k).decay,2)])
    hold off
    xlabel('trial')
    ylabel('slope')
    title(sprintf('Tad %d (k=%d) slope spont trials', allData(k).expnum, k))
    fig_filename = sprintf('Tad %d (k=%d) slope spont trials', allData(k).expnum, k)
    saveas(gcf,fig_filename,'png');
    saveas(gcf, fig_filename, 'epsc2');
    close;
    clear('fig_filename')
end

%%%%%%%%%%%%% no

% Ok, what is the decay constant? And should it be distinct for each tad or
% same across tads?
plot_colors = colormap(spring(length(allData)))
figure;
hold on
for k = 1:length(allData)
    X = k*ones(size(allData(k).decay,2));
    Y = mean(allData(k).decay(:,:,1), 1);
%     if allData(k).stage == 46
%         col = [0 1 0]
%     elseif allData(k).stage == 49
%         col = [1 0 1]
%     else
%         col = [0.5 0.5 0.5]
%     end
    
    plot(X, Y, '*', 'Color', plot_colors(k,:)) %, 'Color', col)
    plot(k, mean(Y), '+', 'MarkerSize', 10, 'Color', [0 0 0], 'LineWidth', 1)
    clear('X', 'Y')
end
hold off
xlim([0 [length(allData) + 1]])
xlabel('Tadpole')
ylabel('Mean slope for each spont trial')
title('Mean slope of each spont trial all tads')
set(gca,'FontSize',20)
fig_filename = 'Mean slope of each spont trial all tads'
saveas(gcf,fig_filename,'png');
saveas(gcf, fig_filename, 'epsc2');

% There is a decent amount of variation across tads, so best to use the grand average for each tad. 

% calculate the function to add to each trace

%%%%%%%%%%%%% leave out for now - concerned about validity of doing this 

%% MISSING FROM IMPORT - ADD NOW

% stimulus onset
% exps 1-9 have onset at 0.5s 
% exps 10+ have onset at 2s
for t = 1:length(allData)
    if allData(t).expnum <=9
        allData(t).stim_onset = 0.5;
    elseif allData(t).expnum >9
        allData(t).stim_onset = 2;
    else
        fprintf('error exp %d, t')
    end
end

% stage
stage(:,1) = [2;3;5;6;7;8;9;10;11;19;20;22;23;24;30;31;32;33;34;35;36;38;40;42;43;44;46;47;48;49];
stage(:,2) = [49;49;49;49;49;49;49;49;49;49;49;49;49;49;46;49;46;1;49;46;46;49;49;49;46;46;46;46;46;46];

for k = 1:length(allData)
    id = find(stage(:,1) == allData(k).expnum)
    allData(k).stage = stage(id,2)
end

%% How many responses? ** Implement stricter criteria here b/c more noise

% calculate response with same 0.1df threshold by trial, but a respROI has
% at least 3 responses in trials that are not spont (4)


for k = 1:length(allData)
    STARTFRAME = floor(allData(k).stim_onset * (allData(k).trial_length(1) / 7)); %stim onset in frames
    ENDFRAME = 2; %how many frames to exclude at end of trial
    allData(k).area_df_f0 = calc_area(allData(k).df_f0, STARTFRAME);
    [allData(k).peak_df_f0, allData(k).peakloc_df_f0] = calc_peak2(allData(k).df_f0, STARTFRAME, ENDFRAME);
    [allData(k).boolean_resp_df_f0, allData(k).sum_resp_all_df_f0] = get_respondingROIs3(allData(k).area_df_f0, allData(k).peak_df_f0, allData(k).peakloc_df_f0);
end

% plot the proportion responses for each tad
figure;
hold on
for k = 1:length(allData)
    X = k * ones(length(allData(k).sum_resp_all_df_f0));
    Y = (allData(k).sum_resp_all_df_f0 ./ length(allData(k).stimorder))
    if allData(k).stage == 46
        col = 'g'
    elseif allData(k).stage == 49
        col = 'm'
    else
        col = 'k'
    end
    plot(X, Y, '*', 'Color', col, 'MarkerSize', 10)
end
hold off
xlim([0 [length(allData) + 1]])
xlabel('Tadpole')
ylabel('Proportion responses')
title('Proportion responses per ROI all tads')
set(gca,'FontSize',20)
fig_filename = 'Proportion responses per ROI all tads'
saveas(gcf,fig_filename,'png');
saveas(gcf, fig_filename, 'epsc2');
    
% proportion non-responsoive ROIs (less than 4 responses)
for k = 1:length(allData)
    non_responders(k) = length(find(allData(k).sum_resp_all_df_f0 < 5)) / length(allData(k).sum_resp_all_df_f0)
end
%Y_vals = 1:1:length(allData);
figure;
bar(non_responders)
xlabel('Tadpole')
ylabel('Proportion unresponsive ROIs')
title('Proportion unresponsive ROIs all tads')
set(gca,'FontSize',20)
fig_filename = 'Proportion unresponsive ROIs all tads'
saveas(gcf,fig_filename,'png');
saveas(gcf, fig_filename, 'epsc2');

% proportion responsive ROIs - at least 1 response that's not in a spont
% trial (stimorder = 4)
for k = 1:length(allData)
    tmp = find(allData(k).stimorder ~= 4);
    tmp1 = sum(allData(k).boolean_resp_df_f0(:,tmp), 2);      
    responders(k) = length(find(tmp1)) / length(tmp1);
    clear('tmp', 'tmp1')
end

% note: at this point, I discovered that k=16 (exp 42) had the wrong number
% of stimuli listed in stimorder. length(stimorder) = 72, but num_trials =
% 60. Removed incorrectly duplicated stimorder infor from 61-72.

%% ID how many good experiments 
% based on 25% or more ROIs having at least 1 with-stimulus trial with a
% peak > 0.1 and area > 0 in the raw calculated df_f0. 

resp_tads = find(responders >= 0.25)
tmp = [allData(:).stage]
ct = tmp(resp_tads)
length(find(ct == 46))
length(find(ct == 49))

% 9 stage 46, 8 stage 49, and the discoboxed one have >= 25% responding
% ROIs

% how many ROIs are in each tad in this group? 
% e.g. is 25% enough to actually do something with?
for k = 1:length(resp_tads)
    tmp = find(allData(resp_tads(k)).stimorder ~= 4);
    tmp1 = sum(allData(resp_tads(k)).boolean_resp_df_f0(:,tmp), 2);      
    respTad_respROIs(k) = length(find(tmp1));
    clear('tmp', 'tmp1') 
end

% yields 11 or more ROIs per exp. So that means we have "enough" ROIs to
% run correlation analysis, etc. 
%%%%%% BUT - N=8 and N=9 is probably too low for publishing purposes. And
%%%%%% some of these ROIs only have 1 response. 

%% Assess and compare ROIs from these 17 "good" experiments to disseration analysis of those tads
% resp_tads has the k value in allData for each of those tads, not the actual exp number. 

%histogram of # of responses for each tad
for k = 1:length(resp_tads)
    tmp = find(allData(resp_tads(k)).stimorder ~= 4);
    tmp1 = sum(allData(resp_tads(k)).boolean_resp_df_f0(:,tmp), 2);      
    [allData(resp_tads(k)).resp_ROIs(:,1), tmp2, allData(resp_tads(k)).resp_ROIs(:,2)] = find(tmp1);
    clear('tmp', 'tmp1') 
end

for k = 1:length(resp_tads)
    total = length(find(allData(resp_tads(k)).stimorder ~= 4));
    num_rois = length(allData(resp_tads(k)).sum_resp_all_df_f0)
    num_resp = length(allData(resp_tads(k)).resp_ROIs)
    num_nonresp = num_rois - num_resp
    figure;
    histogram(allData(resp_tads(k)).resp_ROIs(:,2), ceil(total/3))
    xlabel(sprintf('number of responses (of %d stim pres)', total))
    ylabel('ROI count')
    xlim([0 total])
    %ylim([0 num_resp])
    title(sprintf('Hist responses per respROI tad %d (k=%d)', allData(resp_tads(k)).expnum, resp_tads(k)))
    annotation('textbox', [0.5 0.8 0.1 0.1], 'String', sprintf('%d no-resp ROIs of %d total ROIs', num_nonresp, num_rois'))
    set(gca,'FontSize',20)
    fig_filename = sprintf('Hist responses per respROI tad %d (k=%d)', allData(resp_tads(k)).expnum, resp_tads(k))
    saveas(gcf,fig_filename,'png');
    saveas(gcf, fig_filename, 'epsc2');
    close;
end

%% differences between this analysis and old analysis using smoothed data 

% What ROIs are included? Are they the same ROIs?
% renamed allData to alldata_new and opened AllData_20180320 to be able to
% make the comparison. 
% D:\Torrey_calcium_imaging\compare_46-49\analysis_Feb 2018\corrected_for_badtrials\46_49_comparison

% collect old respROIs into allData_new
for t = 1:length(allData)
    idx_old(t) = allData{1,t}.expnum;
end
idx_new = [allData_new(:).expnum]

for k = 20:length(allData_new)
    t = find(idx_old == idx_new(k))
    allData_new(k).old_respROIs = allData{1,t}.resp_ROIs;
    allData_new(k).ROICenters = allData{1,t}.ROIcenters;
end
% k=9 and k=19 = empty matrix. this corresponds to exp 23 and exp 46 missing from
% allData_new and allData_old contains exps 21 and 28 instead. 21 and 28 do
% not meet criteria for good exp. 
% exp 23 has no respROIs in allData_new. 

% saved allData_new into
% D:\Torrey_calcium_imaging\2018August\allData_versionAug2018_5 and renamed
% to allData when opened.

% How many resp_ROIs old are in resp_ROIs new?
for k = 1:length(allData)
    %if ~isempty(allData(k).resp_ROIs) && ~isempty(allData(k).old_respROIs)
    ct_rois(k,1) = length(setdiff(allData(k).resp_ROIs, allData(k).old_respROIs)); %number rois in new but not old
    ct_rois(k,2) = length(setdiff(allData(k).old_respROIs, allData(k).resp_ROIs)); %number rois in old but not new
    ct_rois(k,3) = length(intersect(allData(k).resp_ROIs, allData(k).old_respROIs)); %number rois in both
    %else
     %   ct_rois(k,:) = [NaN NaN NaN]
    %end
end



%% Summary at end of August 2018:

% Problem 1 - smoothing algorithm created spurious values in the first and
% last frames of a trial. 
% this would increase xcorr b/c I put the trials together. 
% solution is to just go back to using the raw df/f. 

% Problem 2 - some df/f traces had spurious values b/c f0 ~= 0. 
% no solution yet, but could decrease neuropil subtraction. 

% Problem 3 - many ROIs with no activity.
% Activity is much lower than I thought it was, but some tads are ok.
% (resp_tads hold k for allData of these) 

% Dataset info
% allData is a struct with a field for each datapoint.
% The data for allData was imported from: 




