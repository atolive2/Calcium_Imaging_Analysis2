% New and improved calcium imaging analysis script. Uses functions and is
% cool like that.

%% Register Images using imageAnalysis_gui-master from Chris Deister (Moore Lab)
% https://github.com/cdeister/imageAnalysis_gui,
% christopher_deister@brown.edu

% First turn .nd2 files into .tif files and concatenate a trial block
% Unfortunately this must be done 1 by 1, and close files after saving
% them. Change the title to reflect the correct experiment. 

%START BY TRUNCATING ALL TRIALS. Each trial in a block should have the same
%number of frames. 1-2 frames may need to be removed. 
    %Image --> stacks --> tools --> slice remover. "first"=first slice to
    %remove; "last" = last slice in stack, "increment"=1 (so you remove all
    %extra slices)
    %MIJ.run('Slice Remover', 'first=220 last=221 increment=1');
    %selectWindow('seq2610.nd2');
%MIJ.run('Concatenate...', 'all_open title=exp3block12');

% Then Use importer to import each file (multitiff needs to be checked).
% Change the name of importedStack to reflect the trial number

% Then concatenate all trial blocks with the same field of view
flatStack=cat(3, importedStack2, importedStack3);
flatStack=cat(3, flatStack, importedStack3); %repeat for all importedStacks
flatStack=cat(3, flatStack, importedStack4);
flatStack=cat(3, flatStack, importedStack5);
flatStack=cat(3, flatStack, importedStack6);
flatStack=cat(3, flatStack, importedStack7);
flatStack=cat(3, flatStack, importedStack8);
%resulting file should be 500x502xsum(frames) and type uint16 if using 2x2 bin on original
%data

%Then use the GUI to register the images
%Press refresh in the Project/Register column
% select flatStack and then press "luminance" in the Stack Manipulations
% column. This generates a mean luminance profile over all frames. Choose
% about 100 frames near the middle with little to no activity and generate
% a "constrained" projection. If this projection is clear and cells are
% obvious, continue. If not, try a different range of frames

% Once you have a good template, set it as the template (select it, then
% press "set template"). Select flatStack and set it as the registration
% stack (select it, then press "Set Reg St"). Then press "register st"
% which registers th stack of images into flatStack_registered. Generate a
% mean projection to use to make ROIs

% To generate ROIs:
% start roiMaker. Open the meanProj_flatStack_registered by selecting it
% and pressing "Load proje..." Then Press "soma (S o..." to begin ROI
% selection. Zoom in and draw a circle around the ROI, and press "soma (S
% o..." again to add it. If you make a mistake, delete the ROI manually and try again to elimnate NAN ROIs.
%I can make a PCA of the flatStack_registered to
% generate putative ROIs and/or the Global X to make a correlation matrix to help
% select ROIs. The last ROI selected must be a "background" ROI. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Neuropil Extraction is also very important. 
%In roiMaker, enter a value in
%the box labelled "Neuropil PX." This value is the number of pixels that
%the neuropil erodes out from the soma to make a donut shaped area defined
%as the neuropil. It also removes any defined somas from the donut.General
%idea is to have about the same number of pixels in the neuropil as in the
%soma. [current best practice is 18px for a 500x502 image - TT20160726]

% To extract mean fluorescence values across all imags for all ROIs, open roiInspector. 
%Select "somatic" and flatStack_registered. Then press "extract." (not "disk extract") Display
%the mean flourescence over time by selecting ROI Types to Display
%"somatic." This will take a long time - e.g. 7.6mins for 157 ROIs over 12748
%frames (and 20 mins if images are from 2x2 binning). 

%This process leaves you with mean fluoresence for each ROI in each frame
%in somaticF, which holds ROIs in rows and frames in columns (e.g.
%157x12748 double). In addition, you have the soma locations in
%SomaticROICenters. SomaticROICenters is a 1xROI_number cell array. In each
%cell is the X and Y coordinates, stored in Centroid. 0 is in the upper
%left hand corner and it is based on the pixels (so my images will be
%500x502pxls). 

%% Commit basic info to tadpole struct.

% Set up all necessary manual entry variables and values

% 1. stimulus order
    %1 = multisensory high M / high V
    %2 = visual high/crash
    %3 = mechanosensory high
    %4 = no stimulus
    
    %5 = multisensory high M / low V
    %6 = visual low/scrambled crash
    %7 = mechanosensory low
    
    %8 = multisensory low M / high V
    %9 = multisensory low M / low V 
    
    %10 = multisensory med M / high V
    %11 = multisensory med M / low V
    %12 = mechanosensoy med
    
   

% format is 1 long row of length num_trials * numtrialblocks.
% latin square C1 is 1 2 3 4 2 3 1 4 3 1 2 4 (high M / high V)
% latin square C1 is 5 6 3 4 6 3 5 4 3 5 6 4 (high M / low V)
% latin square C1 is 8 2 7 4 2 7 8 4 7 8 2 4 (low M / high V)
% latin square C1 is 9 6 7 4 6 7 9 4 7 5 9 4 (low M / low V)
% latin square C1 is 10 2 12 4 2 12 10 4 12 10 2 4 (med M / high V)
% latin square C1 is 11 6 12 4 6 12 11 4 12 11 6 4 (med M / low V)

% latin square A5 is 3 2 1 4 1 3 2 4 2 1 3 4 (high M /high V)
% latin square A5 is 3 6 5 4 5 3 6 4 6 5 3 4 (high M / low V)
% latin square A5 is 7 2 8 4 8 7 2 4 2 8 7 4 (low M / high V)
% latin square A5 is 7 6 9 4 9 7 6 4 6 9 7 4 (low M / low V)
% latin square A5 is 12 2 10 4 10 12 2 4 2 10 12 4 (med M /high V)
% latin square A5 is 12 6 11 4 11 12 6 4 6 11 12 4 (med M / low V)

% OFFSET EXPS: 
%LS C1 high/high 21 22 23 24 22 23 21 24 23 21 22 24
%LS A5 high/high 23 22 21 24 21 23 22 24 22 21 23 24

% TESTING SHORT VIS (Flash)
% MS = 31; V = 32 M = 3 (same) N = 4 (same)
% LS C1 short(flash)/high 31 32 3 4 32 3 31 4 3 31 32 4
% LS A5 short(flash)/high 3 32 31 4 31 3 32 4 32 31 3 4


% old order is 1 2 3 4 1 2 3 4 1 2 3 4
tadpole.stimorder = [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]
tadpole.stimorder = [1 1 1 1 1 1 1 1 1 1]
tadpole.stimorder = [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]
tadpole.stimorder = [1 2 3 4 2 3 1 4 3 1 2 4 ...
    3 2 1 4 1 3 2 4 2 1 3 4 ...
    9 6 7 4 6 7 9 4 7 5 9 4 ...
    7 6 9 4 9 7 6 4 6 9 7 4 ...
    1 2 3 4 2 3 1 4 3 1 2 4 ...
    3 2 1 4 1 3 2 4 2 1 3 4 ...
    9 6 7 4 6 7 9 4 7 5 9 4 ...
    7 6 9 4 9 7 6 4 6 9 7 4 ...
    1 2 3 4 2 3 1 4 3 1 2 4];
    
    3 2 1 4 1 3 2 4 2 1 3 4 ...
    1 2 3 4 2 3 1 4 3 1 2 4 ...
    3 2 1 4 1 3 2 4 2 1 3 4];
%tadpole.stimorder = [8 2 7 4 2 7 8 4 7 8 2 4 7 2 8 4 8 7 2 4 2 8 7 4 10 2 12 4 2 12 10 4 12 10 2 4 12 2 10 4 10 12 2 4 2 10 12 4 10 2 12 4 2 12 10 4 12 10 2 4]
tadpole.stimorder = [3 2 1 4 1 3 2 4 2 1 3 4] ...
    
7 6 9 4 9 7 6 4 6 9 7 4 ...
    1 2 3 4 2 3 1 4 3 1 2 4 ...
    3 2 1 4 1 3 2 4 2 1 3 4 ...
    5 6 3 4 6 3 5 4 3 5 6 4 ...
    5 6 3 4 6 3 5 4 3 5 6 4]

% alexis stimorder
% 1 = V only 2= M only 3= M-V simul 4=M250V 5=M500V 6 = no stim
[1 2 3 4 5 6 1 2 3 4 5 6];%C1
[6 5 4 3 2 1 6 5 4 3 2 1]; %A1
tadpole.stimorder = [6 5 4 3 2 1 6 5 4 3 2 1 1 2 3 4 5 6 1 2 3 4 5 6];


%check
length(tadpole.stimorder)
%2. experiment number
tadpole.expnum= 0
%3. date of experiment
tadpole.expdate='20180725'
%4. file path
tadpole.filepath= 'L:\\072518\'
%5. make a folder for the figures
mkdir([tadpole.filepath 'figures']); 
%6. trial length
tadpole.trial_length= [160];
%7. number of trial blocks
tadpole.numtrialblocks=2
%8. Create figure save path
tadpole.figure_filepath=[tadpole.filepath 'figures\']
%9. number of trials in a block
tadpole.num_trials=12
%10. what blocks is this?
%tadpole.blockids = [ 1 2 3 4 5 6 7 8 9 ]
%what cell number?
%tadpole.tadid = 1
%tadpole.cellid = 7
%tadpole.cellROI = 1

% %Make sure there are no NaN ROIs--replace all NaN with 0.
% [row, col] = find(isnan(somaticF))
%Manually delete any offending row and record in Analysis notes

%Commit ROI info for somas and neuropil to tadpole.
tadpole.somaticF=somaticF;
%tadpole.somaticF=[somaticF1 somaticF];
tadpole.somaticROI_PixelLists=somaticROI_PixelLists;
tadpole.somaticROIBoundaries=somaticROIBoundaries;
tadpole.somaticROICenters=somaticROICenters;
tadpole.somaticRoiCounter=somaticRoiCounter;
tadpole.somaticROIs=somaticROIs;
tadpole.neuropilF=neuropilF;
%tadpole.neuropilF=[neuropilF1 neuropilF];
tadpole.neuropilROI_PixelLists=neuropilROI_PixelLists;
tadpole.neuropilROIBoundaries=neuropilROIBoundaries;
tadpole.neuropilROICenters=neuropilROICenters;
tadpole.neuropilRoiCounter=neuropilRoiCounter;
tadpole.neuropilROIs=neuropilROIs;

%% At this point, you should save tadpole as a separate .mat file with the filename exp[num]blocks[num-num]_yymmdd.
%save('', '-v7.3')
% then open the mat file containing only tadpole struct.\

% remove extra frames to get 159 frames per trial
% only necessary if forgotten during tiff creation step in image J
% what trials need frames removed?
remove_from_trials = [trial_ids(1:154), trial_ids(156:165), trial_ids(169:175), trial_ids(177), trial_ids(180:186), trial_ids(192:196), trial_ids(198)]

% create a vector to be trial_length
all_trial_lengths = 159 * ones(1,204)
f159_trs = setdiff(trial_ids, remove_from_trials)
for i = 1:length(remove_from_trials)
    all_trial_lengths(remove_from_trials(i)) = all_trial_lengths(remove_from_trials(i)) + 1
end

% split the data into trials
data = tadpole.neuropilF; % switch for somatic and neuropil

frame_end = 0;
frame_start = 1;
num_trials = length(all_trial_lengths)
split_trials=cell(size(data,1),num_trials);
for i = 1:num_trials
    frame_end = frame_end + all_trial_lengths(i)
    for j = 1:size(data,1)
    split_trials{j,i} = data(j, frame_start:frame_end);
    end
    frame_start = frame_end + 1
end
tadpole.trial_splitN_all = split_trials; %switch for somatic and neuropil

% remove last frame of all trials with 160
for i = 1:size(tadpole.trial_splitN_all, 2)
    if length(tadpole.trial_splitN_all{1, i}) > 159
        for j = 1:size(tadpole.trial_splitN_all, 1)
        tadpole.trial_splitN{j, i} = tadpole.trial_splitN_all{j, i}(:, 1:159);
        end
    else
        for j = 1:size(tadpole.trial_splitN_all, 1)
        tadpole.trial_splitN{j, i} = tadpole.trial_splitN_all{j, i};
        end
    end
end

%% Split the data into trials

% start by breaking up trials
[tadpole.trial_splitS] = split_into_trials( tadpole.somaticF, tadpole.trial_length )
[tadpole.trial_splitN] = split_into_trials( tadpole.neuropilF, tadpole.trial_length )

%% Check step: see if there are any terrible trials that should be
% eliminated.
for i = 1:size(tadpole.trial_splitS,2)
    figure;
    hold on
    for j = 1:size(tadpole.trial_splitS,1)
        plot(tadpole.trial_splitS{j,i})
    end
    ax=gca;
    xsize = length(tadpole.trial_splitS{j,i});
    ax.XTick = [0, xsize/7, (xsize/7)*2, (xsize/7)*3, (xsize/7)*4, (xsize/7)*5, (xsize/7)*6, (xsize/7)*7];
    ax.XTickLabel = {'0','1', '2', '3', '4', '5', '6', '7'};
    title(sprintf('Exp %d trial %d all ROIs', tadpole.expnum, i));
    xlabel('time(s)');
    ylabel('raw pixel intensity');
    hold off
    %fig_filename=sprintf([tadpole.figure_filepath, 'exp%dtrial%d.png'], tadpole.expnum, i);
    fig_filename = sprintf('exp%dtrial%d.png', tadpole.expnum, i);
    saveas(gcf,fig_filename,'png');
    close;
    clear('fig_filename')
end

%tadpole.badtrials = [85 86 94 95];

%% After checking for bad trials visually

%tadpole.badtrials = [25, 26, 27, 28];

%% Manipulate "good" trials: 

% Turn raw signal into deltaF/F0
% get a background value for each
[ tadpole.background ] = calc_background( tadpole.trial_splitS )

% subtract background from somaticF and neuropilF
[ tadpole.backgroundsubS ] = subtract_background(tadpole.trial_splitS, tadpole.background)
[ tadpole.backgroundsubN ] = subtract_background(tadpole.trial_splitN, tadpole.background)

% calculate fluorescence decay from somaticF 
%SKIP FOR NOW - MIGHT BE UNNECESSARY

% subtract neuropil from soma to get signal
[ tadpole.signal ] = subtract_neuropil_from_soma( tadpole.backgroundsubS, tadpole.backgroundsubN )

% calculate deltaF/F0
[ tadpole.df_f0 ] = calc_df_f0( tadpole.signal, 10 )

% Extract parameters for each trial
[ tadpole.area_bytrial ] = calc_area( tadpole.df_f0, 42 )
[ tadpole.meanpeak_bytrial, tadpole.peakloc_bytrial] = calc_peak2( tadpole.signal, floor((tadpole.trial_length/7)*2) , 2)
% for unknown reasons, calc_peak suddenly stopped working on gcamp_ours
% (exp 921)
%[ tadpole.meanpeak_bytrial, tadpole.peakloc_bytrial] = calc_peak( tadpole.signal )

% Define response/no response
[ tadpole.boolean_response, tadpole.sum_responses ] = get_respondingROIs( tadpole.area_bytrial, tadpole.meanpeak_bytrial, tadpole.peakloc_bytrial )

% define stim mask by stim type
[ tadpole.stimmask ] = get_stimmask( tadpole.stimorder );
tadpole.unique_stims = unique(tadpole.stimorder);

% Find average area, peak and peakloc for each ROI for each stim type
[ tadpole.area_avg ] = mean_by_stimtype ( cell2mat(tadpole.area_bytrial), tadpole.stimmask )
[ tadpole.peak_avg ] = mean_by_stimtype ( tadpole.meanpeak_bytrial, tadpole.stimmask )
[ tadpole.peakloc_avg ] = mean_by_stimtype ( tadpole.peakloc_bytrial, tadpole.stimmask )

% calculate MS enhancement
[ tadpole.MSenh_area ] = calc_MSenhancement( tadpole.area_avg )
[ tadpole.MSenh_peak ] = calc_MSenhancement( tadpole.peak_avg )
[ tadpole.MSenh_peakloc ] = calc_MSenhancement( tadpole.peakloc_avg )

% collect all values for a single stim type into a single place for easy
% histogram - ing.
[ tadpole.stim_vals_area ] = get_all_bystim( tadpole.stimmask, unique(tadpole.stimorder), cell2mat(tadpole.area_bytrial) )
[ tadpole.stim_vals_meanpeak ] = get_all_bystim( tadpole.stimmask, unique(tadpole.stimorder), tadpole.meanpeak_bytrial )
[ tadpole.stim_vals_peakloc ] = get_all_bystim( tadpole.stimmask, unique(tadpole.stimorder), tadpole.peakloc_bytrial )

%% Smoothed analysis

    for i = 1:size(tadpole.df_f0, 1)
        for j = 1:size(tadpole.df_f0, 2)
            tadpole.smoothed{i,j} = smooth(tadpole.df_f0{i,j}(:,:), 8, 'moving');
        end
    end

%% Recalculate all extracted values using smoothed data


    % Extract parameters for each trial
    [ tadpole.area_bytrial_sm ] = calc_area( tadpole.smoothed, 46 )
    [ tadpole.meanpeak_bytrial_sm, tadpole.peakloc_bytrial_sm] = calc_peak2( tadpole.smoothed, 5, 5)

    % Define response/no response
    [ tadpole.boolean_response_sm, tadpole.sum_responses_sm ] = get_respondingROIs2( tadpole.area_bytrial_sm, tadpole.meanpeak_bytrial_sm, tadpole.peakloc_bytrial_sm )

    % Find average area, peak and peakloc for each ROI for each stim type
    [ tadpole.area_avg_sm ] = mean_by_stimtype2 ( cell2mat(tadpole.area_bytrial_sm), tadpole.stimmask )
    [ tadpole.peak_avg_sm ] = mean_by_stimtype2 ( tadpole.meanpeak_bytrial_sm, tadpole.stimmask )
    [ tadpole.peakloc_avg_sm ] = mean_by_stimtype2 ( tadpole.peakloc_bytrial_sm, tadpole.stimmask )

    % calculate MS enhancement
    % this is only on high/high presentations (1,2,3) because argh
    % calc_MSehancement
    [ tadpole.MSenh_area_sm ] = calc_MSenhancement( tadpole.area_avg_sm )
    [ tadpole.MSenh_peak_sm ] = calc_MSenhancement( tadpole.peak_avg_sm )
    [ tadpole.MSenh_peakloc_sm ] = calc_MSenhancement( tadpole.peakloc_avg_sm )

    % collect all values for a single stim type into a single place for easy
    % histogram - ing.
    [ tadpole.stim_vals_area_sm ] = get_all_bystim( tadpole.stimmask, unique(tadpole.stimorder), tadpole.area_bytrial_sm )
    [ tadpole.stim_vals_meanpeak_sm ] = get_all_bystim( tadpole.stimmask, unique(tadpole.stimorder), tadpole.meanpeak_bytrial_sm )
    [ tadpole.stim_vals_peakloc_sm ] = get_all_bystim( tadpole.stimmask, unique(tadpole.stimorder), tadpole.peakloc_bytrial_sm )
