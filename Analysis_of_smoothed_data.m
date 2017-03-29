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

%% Analysis of peak variability
% calculate difference in peak latency for multisensory vs unisensory
% calculate average peak for each category
