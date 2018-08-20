%% Create master data file with all good tads.

\\files.brown.edu\Research\BM_AizenmanLab\Torrey_calciumimaging

% df/f0, metadata
% put into cell array of structs called raw data

% Open PIE paper assembled data tadpole_all. 

%% Import tads 2-20
for t = 1:length(tads)
    rawData{1,t}.expnum = tads{1,t}.expnum;
    rawData{1,t}.stimorder = tads{1,t}.stimorder;
    rawData{1,t}.df_f0 = tads{1,t}.df_f0;
    for r = 1:length(tads{1,t}.somaticROICenters)
        rawData{1,t}.ROIlocs(r,:) = tads{1,t}.somaticROICenters{1,r}(1).Centroid;
    end
end

%% right data?
t = 1
figure;
hold on
for r = 1:10 %size(tadpole{1,t}.df_f0,1)
    plot(tads.df_f0{r,1})
end
hold off 

%% Find the other tadpoles

% import the ones that are saved as just tadpole (not the whole image extraction workspace)
myFolder = '\\files.brown.edu\Research\BM_AizenmanLab\Torrey_calciumimaging\_Data Analysis\bytad_compile\all_final_tadpolewksp'; % May need to correct this.
if ~isdir(myFolder)
	errorMessage = sprintf('Error: The following folder does not exist:\n%s', myFolder);
	uiwait(warndlg(errorMessage));
	return;
end
filePattern = fullfile(myFolder, 'tadpole*.mat');
matFiles = dir(filePattern)

counter = length(rawData);
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
    %for t = 1:length(matData.tad)
        rawData{1,k+counter} = matData.tadpole % If you get to here, tadpole existed in the file.
        %counter = counter + 1
    %end
    %counter = length(rawData)
end

tads = rawData(1,12:end);
% tads contains exps 21, 22, 23, 24, 30, 31, 32, 34, 40, 42, 43, 45, 46, 48, 49
% some need to have df_f0 extracted from somaticF

%% Find the remaining tadpoles
% 32, 35, 36, 38, 44, 47

% after loading the workspace manually
tads{1, (end + 1)} = tads;

%% What tads do I have?

for t = 1:length(rawData)
    tadIDs(t,1) = rawData{t}.expnum;
end

i = length(tadIDs)
for t = 1:length(tads)
    
    tadIDs((i+t),1) = tads{1,t}.expnum
end

[B, I] = sort(tadIDs)

%% Generate df_f0 for tads that don't have it

for t = 8:length(tads)
    
    [tads{1,t}.trial_splitS] = split_into_trials( tads{1,t}.somaticF, tads{1,t}.trial_length )
    [tads{1,t}.trial_splitN] = split_into_trials( tads{1,t}.neuropilF, tads{1,t}.trial_length )
    % get a background value for each
    [ tads{1,t}.background ] = calc_background( tads{1,t}.trial_splitS )

    % subtract background from somaticF and neuropilF
    [ tads{1,t}.backgroundsubS ] = subtract_background(tads{1,t}.trial_splitS, tads{1,t}.background)
    [ tads{1,t}.backgroundsubN ] = subtract_background(tads{1,t}.trial_splitN, tads{1,t}.background)

    % subtract neuropil from soma to get signal
    [ tads{1,t}.signal ] = subtract_neuropil_from_soma( tads{1,t}.backgroundsubS, tads{1,t}.backgroundsubN )

    % calculate deltaF/F0
    [ tads{1,t}.df_f0 ] = calc_df_f0( tads{1,t}.signal, 10 )
end

%% Add data from tads (exps 21 +) to rawData

for t = 1:length(tads)
    rawData{1,t+11}.expnum = tads{1,t}.expnum;
    rawData{1,t+11}.stimorder = tads{1,t}.stimorder;
    rawData{1,t+11}.df_f0 = tads{1,t}.df_f0;
    for r = 1:length(tads{1,t}.somaticROICenters)
        rawData{1,t+11}.ROIlocs(r,:) = tads{1,t}.somaticROICenters{1,r}(1).Centroid;
    end
end

%% sort rawData by expnum
for t = 1:length(rawData)
    tadIDs(t,1) = rawData{t}.expnum;
end


[B, I] = sort(tadIDs)

for i = 1:length(I)
    rawDataS(1,i) = rawData(1,I(i));
end

for t = 1:length(rawDataS)
    tadIDsS(t,1) = rawDataS{1,t}.expnum;
end

%% Add metadata - stage, stimtime, framerate

% metadata(:,1) contains stage
% metadata(:,2) contains stimtime
% add metadata to rawDataS
for t = 1:length(rawDataS)
    rawDataS{1,t}.stage = metadata(t,1)
    rawDataS{1,t}.stimtime_sec = metadata(t,2)
end

for t = 1:length(rawDataS)
    rawDataS{1,t}.triallength_sec = metadata(t,3)
end

% calculate framerate
%trial length is in metadata(:,3)


%% Add raw signal to rawDataS

% open tadpole_all and make 1 big struct
allTads = [tadpole tads_later];
for t = 1:length(rawData)
    tadIDs_old(t,1) = allTads{t}.expnum;
end
sort(tadIDs_old)

for t = 1:length(rawDataS)
    id = rawDataS{1,t}.expnum
    idx = find(tadIDs_old == id)
    rawDataS{1,t}.rawS = allTads{1,idx}.trial_splitS;
    rawDataS{1,t}.rawN = allTads{1,idx}.trial_splitN;
    rawDataS{1,t}.background = allTads{1,idx}.background;
    rawDataS{1,t}.backgroundsubS = allTads{1,idx}.backgroundsubS;
    rawDataS{1,t}.backgroundsubN = allTads{1,idx}.backgroundsubN;
    rawDataS{1,t}.signal = allTads{1,idx}.signal;
end



