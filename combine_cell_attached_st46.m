%% Combine all cell attached data
% Begin by opening older matfile
\\files.brown.edu\research\BM_AizenmanLab\Torrey_calciumimaging\_Ground Truth\st46_groundtruth
load('st46_allcells.mat')
% data is stored in tadpole
% spike times are in spikeTimes
% one unit in tadpole is 1 cell 
% n = 12 (N = 2)

%% plot each cell's traces (raw df/f0)

for t = 1:length(tadpole)
    figure;
    hold on
    for i = 1:size(tadpole{1,t}.df_f0,2)
        plot(tadpole{1,t}.df_f0{1,i})
    end
    hold off
    ax=gca;
    xsize = length(tadpole{1,t}.df_f0{1,1});
    ax.XTick = [0, xsize/7, (xsize/7)*2, (xsize/7)*3, (xsize/7)*4, (xsize/7)*5, (xsize/7)*6, (xsize/7)*7];
    ax.XTickLabel = {'0','1', '2', '3', '4', '5', '6', '7'};
    title(sprintf('Patched cell %d all trials', t));
    xlabel('time(s)');
    ylabel('raw \DeltaF/F_{0}');
    fig_filename=sprintf(['//files.brown.edu/research/BM_AizenmanLab/Torrey_calciumimaging/_Ground Truth/st46_and_49/' 'Alltrials_df_f0_Cell%d.png'], t);
    saveas(gcf,fig_filename,'png');
    close;
    clear('fig_filename')
end

%% Smooth df/f0

for t = 1:length(tadpole)
    for i = 1:size(tadpole{1,t}.df_f0, 1)
        for j = 1:size(tadpole{1,t}.df_f0, 2)
            tadpole{1,t}.filtered{i,j} = smooth(tadpole{1,t}.df_f0{i,j}(:,:), 8, 'moving');
        end
    end
end

%% Plot each cell's traces, smoothed df/f0
for t = 1:length(tadpole)
    figure;
    hold on
    for i = 1:size(tadpole{1,t}.filtered,2)
        plot(tadpole{1,t}.filtered{1,i})
    end
    hold off
    ax=gca;
    xsize = length(tadpole{1,t}.filtered{1,1});
    ax.XTick = [0, xsize/7, (xsize/7)*2, (xsize/7)*3, (xsize/7)*4, (xsize/7)*5, (xsize/7)*6, (xsize/7)*7];
    ax.XTickLabel = {'0','1', '2', '3', '4', '5', '6', '7'};
    title(sprintf('Patched cell %d all trials smoothed', t));
    xlabel('time(s)');
    ylabel('raw \DeltaF/F_{0}');
    fig_filename=sprintf(['//files.brown.edu/research/BM_AizenmanLab/Torrey_calciumimaging/_Ground Truth/st46_and_49/' 'Alltrials_df_f0_Cell%dsmoothed.png'], t);
    saveas(gcf,fig_filename,'png');
    close;
    clear('fig_filename')
end

%% Begin information extraction from smoothed data
% use regular functions, but on filtered
for t = 1:length(tadpole)
    % Extract parameters for each trial
    [ tadpole{1,t}.area_bytrial_filtered ] = calc_area( tadpole{1,t}.filtered, 42 )
    [ tadpole{1,t}.meanpeak_bytrial_filtered, tadpole{1,t}.peakloc_bytrial_filtered ] = calc_peak2( tadpole{1,t}.filtered, floor((tadpole{1,t}.trial_length/7)*2), 2 )

    % Define response/no response
    [ tadpole{1,t}.boolean_response_filtered, tadpole{1,t}.sum_responses_filtered ] = get_respondingROIs( tadpole{1,t}.area_bytrial_filtered, tadpole{1,t}.meanpeak_bytrial_filtered, tadpole{1,t}.peakloc_bytrial_filtered )

    % define stim mask by stim type
    [ tadpole{1,t}.stimmask_filtered ] = ones(length(tadpole{1,t}.spikeTimes));
    %tadpole.unique_stims = unique(tadpole.stimorder);

    % Find average area, peak and peakloc for each ROI for each stim type
    [ tadpole{1,t}.area_avg_filtered ] = mean_by_stimtype ( tadpole{1,t}.area_bytrial_filtered, tadpole{1,t}.stimmask_filtered )
    [ tadpole{1,t}.peak_avg_filtered ] = mean_by_stimtype ( tadpole{1,t}.meanpeak_bytrial_filtered, tadpole{1,t}.stimmask_filtered )
    [ tadpole{1,t}.peakloc_avg_filtered ] = mean_by_stimtype ( tadpole{1,t}.peakloc_bytrial_filtered, tadpole{1,t}.stimmask_filtered )
end

%% convert spikeTimes in seconds to frame numbers
% 160 frames in 7 sec

for t = 1:length(tadpole)
    for i = 1:length(tadpole{1,t}.spikeTimes)
        tadpole{1,t}.spikeFrames{1,i} = floor(tadpole{1,t}.spikeTimes{1,i} * (160/7));
    end
end

%% get spike count for each trial
for t = 1:length(tadpole)
    for i = 1:length(tadpole{1,t}.spikeTimes)
        tadpole{1,t}.spikeCount{1,i} = length(tadpole{1,t}.spikeTimes{1,i});
    end
end


%% Plot peak vs spike count for all data
figure;
hold on
for t = 1:length(tadpole)
%     if t == 12
%         continue
%     else
    lsp = length(cell2mat(tadpole{1,t}.spikeCount));
    lpk = length(tadpole{1,t}.meanpeak_bytrial_filtered(1,:));
    if lsp ~= lpk 
        spikeCt = [cell2mat(tadpole{1,t}.spikeCount) zeros(1,abs((lsp-lpk)))];
    else
        spikeCt = cell2mat(tadpole{1,t}.spikeCount);
    end
    plot(spikeCt, tadpole{1,t}.meanpeak_bytrial_filtered(1,:), '*')
%     end
end
hold off
% find shitty high spike count cell
legend('1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23', '24', '25', '26', '27', '28')
% shitty cell is 12, which has tons of pre-stim activity.



%% Find trials with early spikes (during baseline)
% earlySpike (before frame 42/1.8 sec) has 1, no early spike has 0.
for t=1:length(tadpole)
    for i = 1:length(tadpole{1,t}.spikeFrames)
        if tadpole{1,t}.spikeFrames{1,i} > 42
            tadpole{1,t}.earlySpikes(1,i) = 0;
        elseif isempty(tadpole{1,t}.spikeFrames{1,i})
            tadpole{1,t}.earlySpikes(1,i) = 0;
        else
            tadpole{1,t}.earlySpikes(1,i) = 1;
        end
    end
end

%% Create a new data structure with every trial separated
% fields t, spikeCt, spikeFrames, earlySpikes, smoothed trace
counter = 1
for t = 1:length(tadpole)
    for i = 1:length(tadpole{1,t}.spikeCount)
        CellsAll46(counter).cellnum = t;
        CellsAll46(counter).trialnum = i;
        CellsAll46(counter).earlySpikes = tadpole{1,t}.earlySpikes(1,i);
        CellsAll46(counter).spikeCount = cell2mat(tadpole{1,t}.spikeCount(1,i));
        CellsAll46(counter).spikeFrames = cell2mat(tadpole{1,t}.spikeFrames(1,i));
        CellsAll46(counter).trace = tadpole{1,t}.filtered{1,i};
        CellsAll46(counter).tracker = [tadpole{1,t}.expnum tadpole{1,t}.cellid];
        counter = counter + 1
    end
end

%% Spike triggered average
% What spike counts do I have traces for?
uniquespikecount = unique([CellsAll46.spikeCount])

% how many traces do I have for each spike Count?
for i = 1:length(uniquespikecount)
    numTraces(i) = sum([CellsAll46.spikeCount] == uniquespikecount(i));
end

% Based on numTraces, only have enough data for grouped traces for 0-8
% spikes. 

% plot all traces for a given spike count (not adjusted for spike time)
for i = 1:12
    figure;
    count = i-1
    idxs = find([CellsAll46.spikeCount] == count)
    hold on
    for j = 1:length(idxs)
        plot(CellsAll46(idxs(j)).trace)
    end
    hold off
    ax=gca;
    xsize = length(tadpole{1,t}.df_f0{1,1});
    ax.XTick = [0, xsize/7, (xsize/7)*2, (xsize/7)*3, (xsize/7)*4, (xsize/7)*5, (xsize/7)*6, (xsize/7)*7];
    ax.XTickLabel = {'0','1', '2', '3', '4', '5', '6', '7'};
    title(sprintf('All trials All cells with %d spikes', count));
    xlabel('time(s)');
    ylabel('smoothed \DeltaF/F_{0}');
    fig_filename=sprintf(['//files.brown.edu/research/BM_AizenmanLab/Torrey_calciumimaging/_Ground Truth/st46_and_49/' 'Alltrials_smoothed_df_f0_%dspikes 46.png'], count);
    saveas(gcf,fig_filename,'png');
    %close;
    clear('fig_filename')
end

% mean and standard deviation for 0 spikes
idxs = find([CellsAll46.spikeCount] == 0)
traces_0spikes = []
for i = 1:length(idxs)
    traces_0spikes = [traces_0spikes CellsAll46(idxs(i)).trace(1:159,1)]
end
avg_0spikes = mean(traces_0spikes,2)
stdev_0spikes = std(traces_0spikes')

%plot mean with std for 0 spikes
figure;
hold on
plot(avg_0spikes, 'k')
plot((avg_0spikes' + stdev_0spikes), 'b')
plot((avg_0spikes' - stdev_0spikes), 'r')
hold off
ax=gca;
xsize = length(avg_0spikes);
ax.XTick = [0, xsize/7, (xsize/7)*2, (xsize/7)*3, (xsize/7)*4, (xsize/7)*5, (xsize/7)*6, (xsize/7)*7];
ax.XTickLabel = {'0','1', '2', '3', '4', '5', '6', '7'};
title('Mean and SD of all traces with 0 spikes');
xlabel('time(s)');
ylabel('smoothed \DeltaF/F_{0}');
fig_filename= 'Mean_and_SD_smoothed_df_f0_0spikes_xadj' %sprintf(['F:/Calcium_Imaging_Analysis/cell_attached_files/Spring2017analysis/figures/' 'Mean_and_SD_smoothed_df_f0_%dspikes.png'], 0);
saveas(gcf,fig_filename,'epsc2');
%close;
clear('fig_filename')

%% Plot based on time of first spike

for i = 2:12
   
    count = i-1
    idxs = find([CellsAll46.spikeCount] == count)
    for k = 1:length(idxs)
        if CellsAll46(idxs(k)).spikeFrames(1,1) > 95 || CellsAll46(idxs(k)).spikeFrames(1,1) < 30
            idxs(k) = 0
        end
    end

    for k = 1:length(idxs)
        if idxs(k) ~= 0
            if min(CellsAll46(idxs(k)).trace(:,:)) < -0.2 
                idxs(k) = 0
            end
        end
    end
    
    figure;
    hold on
    for j = 1:length(idxs)
        if idxs(j) ~= 0 
            start = CellsAll46(idxs(j)).spikeFrames(1,1)
            plot(CellsAll46(idxs(j)).trace((start-30):(start+60), 1))
        end
    end
    plot(31, 0, 'r*')
    hold off
    title(sprintf('Spike Triggered average of traces with %d spikes', count));
    ax=gca;
    xsize = 90;
    ax.XTick = [0, xsize/3.94, (xsize/3.94)*2, (xsize/3.94)*3, (xsize/3.94)*4] %, (xsize/7)*5, (xsize/7)*6, (xsize/7)*7];
    ax.XTickLabel = {'0','1', '2', '3', '4'}; %, '5', '6', '7'};
    ylabel('smoothed \DeltaF/F_{0}');
    xlabel('time(s)');
    fig_filename=sprintf(['//files.brown.edu/research/BM_AizenmanLab/Torrey_calciumimaging/_Ground Truth/st46_and_49/' 'Spike_triggered_avg_smoothed_df_f0_%dspikes.png'], i);
    saveas(gcf,fig_filename,'png');
    close;
    clear('fig_filename')
end
% looks ok, but need to find and eliminate problem traces. 
% problem: they all look pretty similar in size.
% useful note: looks like peak of Ca signal coincides with spike time

%% Figures for Carlos/Chris

%find crap traces
idx = find([CellsAll46.spikeCount] == 4);
for i = 1:length(idx)
    subset2spikes(:,i) = CellsAll46(idx(i)).trace(2:159,1);    
end
[ row col ] = find(subset2spikes(1:20,:) > 0.1)

% for 2 spikes, eliminated 1 trace - exp 18 cell 2, trace 2 (CellsAll(86) )
% for 3 spikes, eliminated 1 trace - exp 27 cell 1, trace 4 (CellsAll(209)

% Average all traces with a given number of spikes
for i = 2:10
    count = i-1;
    idxs = find([CellsAll46.spikeCount] == count);
    for k = 1:length(idxs)
        if CellsAll46(idxs(k)).spikeFrames(1,1) > 95 || CellsAll46(idxs(k)).spikeFrames(1,1) < 30
            idxs(k) = 0;
        end
    end
    idxs_pruned = idxs(idxs~= 0);
    for y = 1:length(idxs_pruned)
        start = CellsAll46(idxs_pruned(y)).spikeFrames(1,1)
        tmpdata(:,y) = CellsAll46(idxs_pruned(y)).trace((start-30):(start+60), 1);
    end
    Mean_traces(i,:) = mean(tmpdata,2);
    clear('idxs', 'idxs_pruned', 'tmpdata')
end
figure; plot(Mean_traces')

% get zero traces
idxs = find([CellsAll46.spikeCount] == 0);
    for y = 1:length(idxs)
        start = 50
        tmpdata(:,y) = CellsAll46(idxs(y)).trace((start-30):(start+60), 1);
    end
Mean_traces(1,:) = mean(tmpdata,2);

% plot average traces (spike triggered averages)
figure; 
hold on
plot(Mean_traces')
plot(31, 0, 'r*')
hold off
title('Spike Triggered averages');
ax=gca;
xsize = 90;
ax.XTick = [0, xsize/3.94, (xsize/3.94)*2, (xsize/3.94)*3, (xsize/3.94)*4] %, (xsize/7)*5, (xsize/7)*6, (xsize/7)*7];
ax.XTickLabel = {'0','1', '2', '3', '4'}; %, '5', '6', '7'};
ylabel('smoothed \DeltaF/F_{0}');
xlabel('time(s)');
legend('0 spikes', '1 spike', '2 spikes', '3 spikes', '4 spikes', '5 spikes', '6 spikes', '7 spike', '8 spikes', 'Spike time', 'Orientation', 'vertical')
fig_filename = 'spike triggered avg 46'
saveas(gcf, fig_filename, 'epsc2')

% plot peak of mean trace by spike count
peaks_spiketrg_avg = max(Mean_traces')
spikecount = [0:8]
%add rational function
spikecount_RF = [0:0.2:8]
p1 = 0.1971
p2 = 0.04119
q1 = 0.8741
for i = 1:length(spikecount_RF)
    rational_fcn(i) = (p1*spikecount_RF(i) + p2)/(spikecount_RF(i) + q1)
end
% rational function conf ints
p1l = 0.1249
p2l = -0.1033
q1l = -1.384
for i = 1:length(spikecount_RF)
    rational_fcnl(i) = (p1l*spikecount_RF(i) + p2l)/(spikecount_RF(i) + q1l)
end

p1u = 0.2693
p2u = 0.1857
q1u = 3.132
for i = 1:length(spikecount_RF)
    rational_fcnu(i) = (p1u*spikecount_RF(i) + p2u)/(spikecount_RF(i) + q1u)
end

figure;
hold on
plot(spikecount, peaks_spiketrg_avg(1:9), 'ok')
%ciplot(rational_fcnl, rational_fcnu, spikecount_RF, 'b'); % this looks
%terrible. 
plot(spikecount_RF, rational_fcn)
hold off
fig_filename = 'peak of mean spike trg avg with rational fcn'
saveas(gcf, fig_filename, 'epsc2')

% add SD of peak to graph
% this looks terrible. 
% 1-8 spikes
for i = 2:9
    count = i-1
    idxs = find([CellsAll46.spikeCount] == count)
    for k = 1:length(idxs)
        if CellsAll46(idxs(k)).spikeFrames(1,1) > 95 || CellsAll46(idxs(k)).spikeFrames(1,1) < 30
            idxs(k) = 0
        end
    end

    for k = 1:length(idxs)
        if idxs(k) ~= 0
            if min(CellsAll46(idxs(k)).trace(:,:)) < -0.2 
                idxs(k) = 0
            end
        end
    end
    for j = 1:length(idxs)
        if idxs(j) ~= 0 
            Peaks_byspikect{i}(j) = max(CellsAll46(idxs(j)).trace(:,:));
        end
    end
end
idxs = find([CellsAll46.spikeCount] == 0);
for y = 1:length(idxs)
    Peaks_byspikect{1}(y) = max(CellsAll46(idxs(y)).trace((start-30):(start+60), 1));
end

for i = 1:length(Peaks_byspikect)
    Peaks_byspikect_avg(i) = mean(Peaks_byspikect{i})
    Peaks_byspikect_sd(i) = std(Peaks_byspikect{i})
end

errorbar(peaks_spiketrg_avg(1:9), Peaks_byspikect_sd(1:9));
ylabel('\DeltaF/F_{0}')
xlabel('spike count')
title('All traces mean with SD')
fig_filename = 'st46 All traces mean with SD by spike ct'
saveas(gcf, fig_filename, 'epsc2')

%% plot spike count vs peak df/f0 for traces with no prestim activity

figure;
spikesvs = [];
dff0vs = [];
hold on 
for c = 1:length(CellsAll46)
    if CellsAll46(c).earlySpikes == 0
        spkct = length(CellsAll46(c).spikeFrames);
        dff0 = max(CellsAll46(c).trace);
        plot(spkct, dff0, 'k*')
        %spikesvs = [spikesvs; spkct];
        %dff0vs = [dff0vs; dff0]
    end
end

% add linear regression
X_pts = 1:10
m = 0.01356
b = 0.1311
plot(X_pts, m*X_pts+b)
hold off
ylabel('\DeltaF/F_{0}')
xlabel('spike count')
title('All traces without pre-stim spikes mean with SD s46')


% get spike count vs dff0 data to use for curve fitting tool
for c = 1:length(CellsAll46)
    if CellsAll46(c).earlySpikes == 0
        spkct_all(c) = length(CellsAll46(c).spikeFrames);
        dff0_all(c) = max(CellsAll46(c).trace);
        %plot(spkct, dff0, 'k*')
        %spikesvs = [spikesvs; spkct];
        %dff0vs = [dff0vs; dff0]
    end
end

% rational function with 1 degree in denom and 1 degree in num looks best. 
% replot with the rational function
X_pts = 0:10
p1 = 0.3545
p2 = 0.2881
q1 = 3.836
for i = 1:length(X_pts)
    rational_fcn(i) = (p1*X_pts(i) + p2)/(X_pts(i) + q1)
end

figure;
hold on
plot(spkct_all, dff0_all, 'k*')
plot(X_pts, rational_fcn)
hold off
ylabel('\DeltaF/F_{0}')
xlabel('spike count')
title('All traces without pre-stim spikes')
fig_filename = 's46 All traces without prestim spikes rational fcn'
saveas(gcf, fig_filename, 'epsc2')

% Plot mean and std dev of each spike count 
% get mean and sd for spike counts with more than 1 trace, if no early
% spikes
clear('tmpdata')
for i = 1:10 %max([CellsAll46.spikeCount])
    count = i-1;
    idxs = find([CellsAll46.spikeCount] == count);
    if length(idxs) > 1
        for j = 1:length(idxs)
            if CellsAll46(idxs(j)).earlySpikes == 0
            tmpdata(:,j) = [CellsAll46(idxs(j)).trace(1:159, 1)];
            end
        end
        if exist('tmpdata', 'var')
        max_tmpdata = max(tmpdata); 
        length(max_tmpdata)
        max_tmpdata(max_tmpdata == 0) = [];
        length(max_tmpdata)
        avg_dff0(i) = mean(max_tmpdata);
        sddev_dff0(i) = std(max_tmpdata);
        clear('tmpdata', 'max_tmpdata');
        end
    end
end

% plot mean and std
figure;
errorbar(avg_dff0, sddev_dff0, 'o')
ax=gca;
xsize = length(avg_dff0);
xlim([0 11])
ax.XTick = [1, 6, 11, 16, 21, 26];;
ax.XTickLabel = {'0','5', '10', '15', '20', '25'}; 
ylabel('mean max \DeltaF/F_{0}')
xlabel('spike count')
title('spike count vs mean max \DeltaF/F_{0}')

%% Example cell plots

% Cell 10
figure;
plot([CellsAll46(52:56).trace])
legend('3', '0', '1', '2', '3')

% cell 9
figure; 
plot([CellsAll46(38:51).trace])
legend('17','2','13','8','0','8','2','0','7','1','3','2','3','1')
spikenums = unique([CellsAll46(38:51).spikeCount])
figure;
hold on
for c = 38:51
    if CellsAll46(c).spikeCount == 0
        plot(CellsAll46(c).trace, 'k')
    elseif CellsAll46(c).spikeCount == 1
        plot(CellsAll46(c).trace, 'b')
    elseif CellsAll46(c).spikeCount == 2
        plot(CellsAll46(c).trace, 'r')  
    elseif CellsAll46(c).spikeCount == 3
        plot(CellsAll46(c).trace, 'm')
    elseif CellsAll46(c).spikeCount > 6
        plot(CellsAll46(c).trace, 'Color', [0 1 1])
    end
end
hold off
    ax=gca;
    xsize = 160;
    ax.XTick = [0, xsize/7, (xsize/7)*2, (xsize/7)*3, (xsize/7)*4, (xsize/7)*5, (xsize/7)*6, (xsize/7)*7];
    ax.XTickLabel = {'0','1', '2', '3', '4', '5', '6', '7'};
    ylabel('smoothed \DeltaF/F_{0}');
    xlabel('time(s)');
for c = 38:51
    toexport = mat2str(CellsAll46(c).spikeFrames)
end

%% Plot all with prestim spikes
figure;
hold on
for c = 1:length(CellsAll46)
    if CellsAll46(c).earlySpikes == 1
        plot(CellsAll46(c).trace)
    end
end
hold off
    ax=gca;
    xsize = 160;
    ax.XTick = [0, xsize/7, (xsize/7)*2, (xsize/7)*3, (xsize/7)*4, (xsize/7)*5, (xsize/7)*6, (xsize/7)*7];
    ax.XTickLabel = {'0','1', '2', '3', '4', '5', '6', '7'};
    ylabel('smoothed \DeltaF/F_{0}');
    xlabel('time(s)');
    title('all with early spikes')
    
% find and remove bad traces
earlytraces = find([CellsAll46(:).earlySpikes] == 1)
[cell loc] = ([CellsAll46(earlytraces).trace] > 0.3)

for c = 1:length(earlytraces)
    tmpdata(:,c) = [CellsAll46(c).trace(1:159, 1)];
end

max(tmpdata)
% shitty trace is 207
min(tmpdata)
% shitty trace is 216

% plot only useful traces
figure;
hold on
for c = 1:length(earlytraces)
        plot(CellsAll46(earlytraces(c)).trace)
        %pause(1)
end

hold off
    ax=gca;
    xsize = 160;
    ax.XTick = [0, xsize/7, (xsize/7)*2, (xsize/7)*3, (xsize/7)*4, (xsize/7)*5, (xsize/7)*6, (xsize/7)*7];
    ax.XTickLabel = {'0','1', '2', '3', '4', '5', '6', '7'};
    ylabel('smoothed \DeltaF/F_{0}');
    xlabel('time(s)');
    title('all with early spikes')