%% Stage 46 ground truth analysis
% N = 3, n = 12
% mostly copied from combine_cell_attached, which was used for stage 49
% ground truth data

%% combine all data into 1 mat file
myFolder = 'D:/Torrey_calcium_imaging/st46_groundtruth/'; % May need to correct this.
mkdir([myFolder 'figures']);
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
end

%% Import spike times 

load('st46_spiketimes.mat')

% import spike times by tadpole
for t = 1:length(tadpole)
    tadpole{1,t}.spikeTimes = spikeTimes{t};
end

%% Split data by trial, turn into df/f0 and smooth it

% start by breaking up trials
for t = 1:length(tadpole)
    [tadpole{1,t}.trial_splitS] = split_into_trials( tadpole{1,t}.somaticF, tadpole{1,t}.trial_length )
    [tadpole{1,t}.trial_splitN] = split_into_trials( tadpole{1,t}.neuropilF, tadpole{1,t}.trial_length )
end

%% Check step: see if there are any terrible trials that should be
% eliminated.
for t = 1:length(tadpole)
    for i = 1:size(tadpole{1,t}.trial_splitS,2)
        figure;
        hold on
        for j = 1:size(tadpole{1,t}.trial_splitS,1)
            plot(tadpole{1,t}.trial_splitS{j,i})
        end
        ax=gca;
        xsize = length(tadpole{1,t}.trial_splitS{j,i});
        ax.XTick = [0, xsize/7, (xsize/7)*2, (xsize/7)*3, (xsize/7)*4, (xsize/7)*5, (xsize/7)*6, (xsize/7)*7];
        ax.XTickLabel = {'0','1', '2', '3', '4', '5', '6', '7'};
        title(sprintf('Cell %d trial %d all ROIs', t, i));
        xlabel('time(s)');
        ylabel('raw pixel intensity');
        hold off
        fig_filename=sprintf(['D:/Torrey_calcium_imaging/st46_groundtruth/figures/cell%dtrial%d.png'],t, i);
        saveas(gcf,fig_filename,'png');
        close;
        clear('fig_filename')
    end
end

%% % Turn raw signal into deltaF/F0
for t = 1:length(tadpole)
    % get a background value for each
    [ tadpole{1,t}.background ] = calc_background( tadpole{1,t}.trial_splitS )

    % subtract background from somaticF and neuropilF
    [ tadpole{1,t}.backgroundsubS ] = subtract_background(tadpole{1,t}.trial_splitS, tadpole{1,t}.background)
    [ tadpole{1,t}.backgroundsubN ] = subtract_background(tadpole{1,t}.trial_splitN, tadpole{1,t}.background)

    % subtract neuropil from soma to get signal
    [ tadpole{1,t}.signal ] = subtract_neuropil_from_soma( tadpole{1,t}.backgroundsubS, tadpole{1,t}.backgroundsubN )

    % calculate deltaF/F0
    [ tadpole{1,t}.df_f0 ] = calc_df_f0( tadpole{1,t}.signal, 10 )

    % smooth deltaF/F0
    for i = 1:size(tadpole{1,t}.df_f0, 1)
        for j = 1:size(tadpole{1,t}.df_f0, 2)
            tadpole{1,t}.smoothed{i,j} = smooth(tadpole{1,t}.df_f0{i,j}(:,:), 8, 'moving');
        end
    end
end

%% extract basic information about trials
for t = 1:length(tadpole)
    [ tadpole{1,t}.area_bytrial_sm ] = calc_area( tadpole{1,t}.smoothed, 46 )
    [ tadpole{1,t}.meanpeak_bytrial_sm, tadpole{1,t}.peakloc_bytrial_sm] = calc_peak2( tadpole{1,t}.smoothed, 5, 5)

    % Define response/no response
    [ tadpole{1,t}.boolean_response_sm, tadpole{1,t}.sum_responses_sm ] = get_respondingROIs2( tadpole{1,t}.area_bytrial_sm, tadpole{1,t}.meanpeak_bytrial_sm, tadpole{1,t}.peakloc_bytrial_sm )
end

%% Plot each of the patched cells smoothed df/f0

% only include trials with spikes after stim only
%t=11
% what is max color range?
all_stimcounts = [];
for t = 1:length(tadpole)
    for i = 1:tadpole{1,t}.num_trials
        if ~isempty(tadpole{1,t}.spikeTimes{i})
            if tadpole{1,t}.spikeTimes{i}(1) > artifact_s
                all_stimcounts = [all_stimcounts, length(tadpole{1,t}.spikeTimes{i})];
             end
        end
    end
end
max_spikect = max(all_stimcounts)
hist(all_stimcounts)

line_colors = jet(6)
for t = 1:length(tadpole)
    figure;
    hold on
    for i = 1:tadpole{1,t}.num_trials
        if ~isempty(tadpole{1,t}.spikeTimes{i})
            if tadpole{1,t}.spikeTimes{i}(1) > artifact_s
                line_color = line_colors(length(tadpole{1,t}.spikeTimes{i}), :)
                plot(tadpole{1,t}.smoothed{tadpole{1,t}.cellROI, i}, 'Color', line_color)
            end
        end
    end
    ax=gca;
    xsize = length(tadpole{1,t}.smoothed{tadpole{1,t}.cellROI, 1});
    ax.XTick = [0, xsize/7, (xsize/7)*2, (xsize/7)*3, (xsize/7)*4, (xsize/7)*5, (xsize/7)*6, (xsize/7)*7];
    ax.XTickLabel = {'0','1', '2', '3', '4', '5', '6', '7'};
    title(sprintf('Cell %d stim-locked response trials only', t));
    xlabel('time(s)');
    ylabel('smoothed \DeltaF/F_{0}');
    hold off
    fig_filename=sprintf(['D:/Torrey_calcium_imaging/st46_groundtruth/figures/stim-locked response trials only cell%d.png'],t);
    saveas(gcf,fig_filename,'png');
    close;
    clear('fig_filename')
end

%% Combine across cells: plot peak vs number of spikes 
peaks_byspikect = cell(1,(max_spikect+1));
for t = 1:length(tadpole)
    for i = 1:length(tadpole{1,t}.spikeTimes)
        if ~isempty(tadpole{1,t}.spikeTimes{i})
            if tadpole{1,t}.spikeTimes{i}(1) > artifact_s
                j = length(tadpole{1,t}.spikeTimes{i})
                    peaks_byspikect{1,(j+1)} = [peaks_byspikect{1, j+1}, tadpole{1,t}.meanpeak_bytrial_sm(tadpole{1,t}.cellROI, i)]
    
            end
        else
           peaks_byspikect{1,1} = [peaks_byspikect{1, 1}, tadpole{1,t}.meanpeak_bytrial_sm(tadpole{1,t}.cellROI, i)]
        end
    end
end

% scatterplot
figure;
hold on
for i = 1:size(peaks_byspikect, 2)
    x_vals = i * ones(1, length(peaks_byspikect{1,i}));
    scatter(x_vals, peaks_byspikect{1,i})
    clear('x_vals')
end
hold off
ax=gca;
ax.XTick = [1:1:length(peaks_byspikect)];
ax.XTickLabel = [0:1:(length(peaks_byspikect)-1)];
title('Peak \DeltaF/F_{0} by spike count')
xlabel('total number of spikes')
ylabel('peak \DeltaF/F_{0}')
fig_filename = 'Peak by total number of spikes st46'
saveas(gcf,fig_filename,'png');
close;
clear('fig_filename')

%% Average each df/f0 trace by spike count and plot

% get all spike counts into a vector
all_spike_nums = [];
for t = 1:length(tadpole)
    for i = 1:tadpole{1,t}.num_trials
        if ~isempty(tadpole{1,t}.spikeTimes{i})
            if tadpole{1,t}.spikeTimes{i}(1) > artifact_s
                all_spike_nums = [all_spike_nums, length(tadpole{1,t}.spikeTimes{i})];
            else
                all_spike_nums = [all_spike_nums, -1];
            end
        else
            all_spike_nums = [all_spike_nums, 0];
        end
    end
end

% get all smoothed df/f0 traces into a matrix
all_dff0 = [];
for t = 1:length(tadpole)
    for i = 1:tadpole{1,t}.num_trials
        all_dff0 = [all_dff0, tadpole{1,t}.smoothed{tadpole{1,t}.cellROI, i}];
    end
end

spike_counts = unique(all_spike_nums)

for i = 1:length(spike_counts)
    trials = all_spike_nums == spike_counts(i)
        sorted_dff0{i} = all_dff0(:,trials)
end

% average by spike count
for i = 1:length(sorted_dff0)
    avg_trials(:,i) = mean(sorted_dff0{i},2);
end

plot_colors = jet(size(avg_trials,2))
figure;
hold on
for i = 1:size(avg_trials,2)
    plot(avg_trials(:,i), 'Color', plot_colors(i,:))
end
hold off
colorbar('TickLabels', spike_counts)
ax=gca;
xsize = length(tadpole{1,t}.smoothed{tadpole{1,t}.cellROI, 1});
ax.XTick = [0, xsize/7, (xsize/7)*2, (xsize/7)*3, (xsize/7)*4, (xsize/7)*5, (xsize/7)*6, (xsize/7)*7];
ax.XTickLabel = {'0','1', '2', '3', '4', '5', '6', '7'};
title(sprintf('Average of all stim-locked response trials by spike count', t));
xlabel('time(s)');
ylabel('smoothed \DeltaF/F_{0}');
fig_filename=sprintf('D:/Torrey_calcium_imaging/st46_groundtruth/figures/Average of all stim-locked response trials by spike count stg46.png');
saveas(gcf,fig_filename,'png');
close;
clear('fig_filename')

