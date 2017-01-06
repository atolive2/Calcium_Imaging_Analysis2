%% Loop through all experiments to generate PIE based on peak_avg 
% with responses < 0.1 in both uni and multi eliminated.
% find files of interest
myFolder = 'F:/Calcium_Imaging_Analysis/tadpoles_byexp/'; % May need to correct this.
mkdir([myFolder 'figures']);
if ~isdir(myFolder)
	errorMessage = sprintf('Error: The following folder does not exist:\n%s', myFolder);
	uiwait(warndlg(errorMessage));
	return;
end
filePattern = fullfile(myFolder, 'exp*.mat');
matFiles = dir(filePattern)

% get all tadpoles into a struct
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
	tadpole{k} = matData.tadpole; % Extract the tadpole from the structure.
end


for t = 1:length(tadpole) 
    % for each tadpole, generate max uni and max multi vals
    [tadpole{1,t}.unimax_peakavg, tadpole{1,t}.unimax_stimtype] = max(tadpole{1,t}.peak_avg(2:3,:),[],1);
    tadpole{1,t}.multimax_peakavg = tadpole{1,t}.peak_avg(1,:);
    
    % for each tadpole, generate responding ROIs ( max uni or multi > 0.1)
    tadpole{1,t}.respROIs = find(tadpole{1,t}.unimax_peakavg > 0.1 | tadpole{1,t}.multimax_peakavg > 0.1);
    sprintf('tadpole %d has %d responding ROIs', t, length(tadpole{1,t}.respROIs))
end

% for each tadpole, calculate the exponential fit
for t = 1:length(tadpole)
    tadpole{1,t}.f1 = fit(tadpole{1,t}.unimax_peakavg(1,tadpole{1,t}.respROIs)', tadpole{1,t}.multimax_peakavg(:, tadpole{1,t}.respROIs)', 'exp1'); % this is a 1 term model
    if length(tadpole{1,t}.respROIs) > 4 % 2 term exp requires at least 4 points
        tadpole{1,t}.f2 = fit(tadpole{1,t}.unimax_peakavg(1,tadpole{1,t}.respROIs)', tadpole{1,t}.multimax_peakavg(:, tadpole{1,t}.respROIs)','exp2'); % this is a 2 term model
    end
end

% for each tadpole, plot max uni vs multi by color with exp fit
%%%%%%% this doesn't work yet - need to incorporate delinating max uni type
for t = 1:length(tadpole)
   uniV = 
   figure;
   plot(
end

% example from old code
[uni_largestval modality] = max(largest_vals(:,2:3),[],2)
multi_largestval_toplot = largest_vals(:,1)
primary_vis =[];
primary_mech = [];
for i = 1:length(modality)
    if modality(i)==1
        primary_vis = [primary_vis; uni_largestval(i) multi_largestval_toplot(i)]
    elseif modality(i)==2
        primary_mech = [primary_mech; uni_largestval(i) multi_largestval_toplot(i)]
    else
        sprintf('error')
    end
end


%% Collapse across all tadpoles

% generate 1 matrix with all respROIs from all tads. 
% 1 = multi peak_avg
% 2 = max uni peak_avg
% 3 = modality of max uni (1 = vis, 2 = mech)
% 4 = tadpole ID (10, 11, 19, 20, 2, 3, 5, 6, 7, 8, 9)

pie_allrespROIs = [];
for t = 1:length(tadpole)
    tmpdata = [tadpole{1,t}.multimax_peakavg(1, tadpole{1,t}.respROIs); tadpole{1,t}.unimax_peakavg(1, tadpole{1,t}.respROIs); tadpole{1,t}.unimax_stimtype(1, tadpole{1,t}.respROIs)]
    tadID = t * ones(1,size(tmpdata,2))
    pie_allrespROIs = [pie_allrespROIs, [tmpdata; tadID]];
end

figure;
plot(pie_allrespROIs(1,:), pie_allrespROIs(2,:), 'o')
axis([-0.1 3 -0.1 3])

% eliminated 3 ROIs - 126, 170 and 131 based on uni or multi > 4

% plot uni max vs MS enhancement (that's what shows PIE in my old ppt)
% for all responding ROIs from all tadpoles, calculate a grand exp fit and
% plot
MSenh = multi ./ uni
fit2 = fit(uni', MSenh', 'exp2');

figure;
hold on
plot(uni, MSenh, 'ko')
plot(fit2)
hold off
axis([-0.1 2.5 -1 40])
xlabel('max unisensory response', 'fontsize', 20)
ylabel('multisensory enhancement', 'fontsize', 20)
set(gca,'FontSize',20)
%title('Inverse effectiveness, all responding ROIs for all tadpoles')

% create inset
figure;
hold on
plot(uni, MSenh, 'ko')
%plot(fit1)
hold off
axis([0 0.5 -1 10])
xlabel('max unisensory response', 'fontsize', 20)
ylabel('multisensory enhancement', 'fontsize', 20)
set(gca,'FontSize',20)

