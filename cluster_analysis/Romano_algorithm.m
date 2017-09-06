%% Implement Romano clustering algorithm

% first open the workpace and get the mean imag
% F:\Calcium_Imaging_Analysis\analyzed_compiled\compile_exp6\field1_wprkspace
save('exp6_meanimg.mat', 'meanProj_flatStack_registered')
%then clear workspace and open the mean image and the tadpole data file

% begin by opening the tadpole data file
%load('F:\Calcium_Imaging_Analysis\analyzed_compiled\Smoothed_analysis\exp6tadpole_blocks1-8.mat')

% generate required datafile
% exp6_RASTER.mat

%mean image (from other file)
dataAllCells.avg = meanProj_flatStack_registered

%deltaFoF = T imaging frames x N ROIs matrix of df/f0
% recombine trial by trial df/f0 into 1 long vector for each ROI

    for i = 1:size(tadpole.df_f0, 1)
        tmp = [];
        for j = 1:size(tadpole.df_f0, 2)
            tmp = [tmp tadpole.df_f0{i,j}];
        end
        deltaFoF(:,i) = tmp';
    end

%raster: a TxN matrix of one's to include all frames
raster = ones(size(deltaFoF));

% movements: TxN matrix with 1's for aritifact containing frames, 0's
% elsewhere
movements = zeros(size(deltaFoF));

% perimeter coordinates of each ROI
% clean the data first - some are 1x1 and some are 2x1. Keep cell{1,1}. 
for i = 1:length(tadpole.somaticROIBoundaries)
    if sum(size(tadpole.somaticROIBoundaries{1,i})) ~= 2
        sprintf('doubles in %d', i)
        tadpole.somaticROIBoundaries{1,i}{2,1} = []; %tadpole.somaticROIBoundaries{1,i}{1,1}
    end
end
dataAllCells.cell_per = tadpole.somaticROIBoundaries';

% pixel indexes of each ROI
dataAllCells.cell = tadpole.somaticROI_PixelLists

%save these items to exp6_RASTER.mat
save('exp6_RASTER.mat', 'deltaFoF', 'raster', 'movements', 'dataAllCells')






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% response analysis input file
% exp6_responseAnalysisInput
%1xS cell array called mapData. S=total number of experimental event types
S = unique(tadpole.stimorder)
mapData = num2cell(S)
% labels
for i = 1:length(mapData)
    mapData{1,i}.label = num2str(S(i))
    mapData{1,i}.value = S(i)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%onset time - finish this later, only if I want response analysis.
%get trial lengths in a vector

if length(tadpole.trial_length) == 1
    trial_length = tadpole.trial_length *ones(1, numtrialblocks)
else
    trial_length = tadpole.trial_length
end
all_trial_lengths = reshape(trial_length(ones(12,1),:), 1, [])

%get stim onset for each trial
for i = 1:length(mapData)
    trials = find(tadpole.stimorder(S(i)))
    for t = 1:length(trials)
        onset_time(t) = 
    end
end

    mapData{1,i}.onsetTime = 

%%%%%%%%%%%%%%%%%%%%%%






dataAllCells.avg