function [ tadpole ] = get_matFiles( myFolder, matFiles )
%get_matFiles retrieves tadpole from exp*.mat files from a given folder.
%   input: matFiles, a list of filenames from a given folder
%   output: data, a 1 by length(matFiles) cell array containing all the data
%   from tadpole in each matFile.

for k = 1:length(matFiles)
	matFilename = fullfile(myFolder, matFiles(k).name);
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
	tadpole{k} = matData.tadpole; % Extract the tadpole from the structure.
end
end

