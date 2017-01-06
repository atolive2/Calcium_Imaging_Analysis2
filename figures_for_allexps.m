%% For each experiment, make figures. 

myFolder = 'F:/Calcium_Imaging_Analysis/tadpoles_byexp/'; % May need to correct this.
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
	% If you get to here, tadpole existed in the file.
	tadpole = matData.tadpole; % Extract the tadpole from the structure.
    [unisensory, uniloc] = max(tadpole.area_avg(2:3,:),[],1);
    figure;
    plot(unisensory,tadpole.MSenh_area(1,:), '*');
    title(sprintf([matFiles(k).name, ' area-based MS enhancement raw values']));
    xlabel('largest unisensory');
    ylabel('multisensory enhancement');
    fig_filename=(sprintf([myFolder, 'figures/' matFiles(k).name]));
    saveas(gcf,fig_filename,'png');
    % 	for row = 1 : size(tadpole,1)
% 		% Calculate some stuff
% 	end
end