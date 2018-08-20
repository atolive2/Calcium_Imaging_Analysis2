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





