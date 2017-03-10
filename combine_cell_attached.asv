%% Combine all cell attached data
%alignment of cells with tadpole #
%1=e12c3
%2=e12c4
%3=e12c6
%4=e13c1
%5=e14c1
%6=e14c2
%7=e14c3
%8=e14c4
%9=e15c1
%10=e15c2
%11=e17c1
%12=e17c2
%13=e18c1
%14=e18c3
%15=e18c4
%16=e25c1
%17=e25c2
%18=e26c1
%19=e26c2
%20=e27t1c1
%21=e27t2c1
%22=e27t2c2

%% combine all data into 1 mat file
myFolder = 'F:/Calcium_Imaging_Analysis/cell_attached_files/Spring2017analysis/bycell_data/'; % May need to correct this.
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

%% import spike times (from earlier analysis)
for t = 1:15
    tadpole{1,t}.spikeTimes = AllCells(t).spikeTimes
end

