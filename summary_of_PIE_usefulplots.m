%% summary of PIE all exps just MSenh and LS

myFolder = 'F:/Calcium_Imaging_Analysis/tadpoles_byexp/'; % May need to correct this.
mkdir([myFolder 'figures']);
if ~isdir(myFolder)
	errorMessage = sprintf('Error: The following folder does not exist:\n%s', myFolder);
	uiwait(warndlg(errorMessage));
	return;
end
filePattern = fullfile(myFolder, 'exp*.mat');
matFiles = dir(filePattern)

area_V =[];
area_M = [];
area_Vls = [];
area_Vm = [];
area_Mls = [];
area_Mm = [];
peak_V = []; 
peak_M = []; 
peak_Vls = [];
peak_Vm = [];
peak_Mls = [];
peak_Mm = [];


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
    %subplot 1 - area: lgst uni vs MS enh
    [unisensory, uniloc] = max(tadpole.area_avg(2:3,:),[],1); %get val and loc of max unisensory
    area_V = horzcat(area_V, [unisensory(uniloc==1); tadpole.area_avg(uniloc==1)]); % where max unisensory = V (1)
    area_M = horzcat(area_M, [unisensory(uniloc==2); tadpole.area_avg(uniloc==2)]); % where max unisensory = M (2)

    % subplot 2: area - linear sum vs MS response
    % get data to plot 
    area_Vls = horzcat(area_Vls, tadpole.area_avg(2,uniloc==1) + tadpole.area_avg(3,uniloc == 1));
    area_Vm = horzcat(area_Vm, tadpole.area_avg(1,uniloc==1));
    area_Mls = horzcat(area_Mls, tadpole.area_avg(2,uniloc==2) + tadpole.area_avg(3,uniloc == 2));
    area_Mm = horzcat(area_Mm, tadpole.area_avg(1,uniloc==2));
      
    % subplot 3: peak - lgst uni vs MSenh all
    % get data to plot, by type of max unisensory
    [unisensory, uniloc] = max(tadpole.peak_avg(2:3,:),[],1); %get val and loc of max unisensory
    peak_V = horzcat(peak_V, [unisensory(uniloc==1); tadpole.peak_avg(uniloc==1)]); % where max unisensory = V (1)
    peak_M = horzcat(peak_M, [unisensory(uniloc==2); tadpole.peak_avg(uniloc==2)]); % where max unisensory = M (2)

    % subplot 4: peak - linear sum vs MS response
    % get data to plot 
    peak_Vls = horzcat(peak_Vls, tadpole.peak_avg(2,uniloc==1) + tadpole.peak_avg(3,uniloc == 1));
    peak_Vm = horzcat(peak_Vm, tadpole.peak_avg(1,uniloc==1));
    peak_Mls = horzcat(peak_Mls, tadpole.peak_avg(2,uniloc==2) + tadpole.peak_avg(3,uniloc == 2));
    peak_Mm = horzcat(peak_Mm, tadpole.peak_avg(1,uniloc==2));

end

figure;
% %subplot 1: area MSenh
% % get values for plotting Y=X line
% uni_mm = [floor(min(min([area_V(1,:) area_M(1,:)]))), ceil(max(max([area_V(1,:), area_M(1,:)])))];
% multi_mm = [floor(min(min([area_V(2,:), area_M(2,:)]))), ceil(max(max([area_V(2,:), area_M(2,:)])))];
% Y = (min([uni_mm, multi_mm])-1):1:(max(([uni_mm, multi_mm]))+1);
% 
% subplot(2,2,1)
% hold on
% plot(area_V(1,:), area_V(2,:), 'go')
% plot(area_M(1,:), area_M(2,:), 'mo')
% plot(Y, Y, 'b-')
% ylim([-10 100])
% xlim([-10 150])
% title('Area: largest uni vs MS enh')
% xlabel('largest unisensory mean')
% ylabel('MS enhancement mean')
% legend('Vis', 'Mech', 'Y=X', 'location', 'best')
% hold off
% 
% % subplot 2: area - linear sum vs MS response
% % get values for plotting Y=X line
% uni_mm_ls = [floor(min(min([area_Vls(1,:) area_Mls(1,:)]))), ceil(max(max([area_Vls(1,:), area_Mls(1,:)])))];
% multi_mm_m = [floor(min(min([area_Vm(1,:), area_Mm(1,:)]))), ceil(max(max([area_Vm(1,:), area_Mm(1,:)])))];
% Y = (min([uni_mm_ls, multi_mm_m])-1):1:(max(([uni_mm_ls, multi_mm_m]))+1);
% 
% subplot(2,2,2)
% hold on
% plot(area_Vls(1,:), area_Vm(1,:), 'go')
% plot(area_Mls(1,:), area_Mm(1,:), 'mo')
% plot(Y, Y, 'b-')
% ylim([-75 100])
% xlim([-75 150])
% title('Area: linear sum vs ms')
% xlabel('linear sum mean')
% ylabel('multisensory mean')
% legend('Vis', 'Mech', 'Y=X', 'location', 'best')
% hold off

%subplot 3: peak - lgst uni vs MSenh all
% get values for plotting Y=X line
uni_mm = [min(min([peak_V(1,:) peak_M(1,:)])), max(max([peak_V(1,:), peak_M(1,:)]))];
multi_mm = [min(min([peak_V(2,:), peak_M(2,:)])), max(max([peak_V(2,:), peak_M(2,:)]))];
Y = (min([uni_mm, multi_mm])):0.01:(max(([uni_mm, multi_mm])));

subplot(1,2,1)
hold on
plot(peak_V(1,:), peak_V(2,:), 'ro', 'MarkerFaceColor', 'r')
plot(peak_M(1,:), peak_M(2,:), 'bo', 'MarkerFaceColor', 'b')
%plot(Y, Y, 'b-')
ylim([-0.01 0.7])
xlim([-0.01 2.2])
title('Largest Unisensory Response vs Multisensory Enhancement')
xlabel('largest unisensory mean')
ylabel('multisensory enhancement mean')
legend('Vis', 'Mech', 'Y=X', 'location', 'best')
hold off

% subplot 4: peak - linear sum vs MS response
% get values for plotting Y=X line
uni_mm_ls = [min(min([peak_Vls(1,:) peak_Mls(1,:)])), max(max([peak_Vls(1,:), peak_Mls(1,:)]))];
multi_mm_m = [min(min([peak_Vm(1,:), peak_Mm(1,:)])), max(max([peak_Vm(1,:), peak_Mm(1,:)]))];
Y = (min([uni_mm_ls, multi_mm_m])):0.01:(max(([uni_mm_ls, multi_mm_m])));

subplot(1,2,2)
hold on
plot(peak_Vls(1,:), peak_Vm(1,:), 'ro', 'MarkerFaceColor', 'r')
plot(peak_Mls(1,:), peak_Mm(1,:), 'bo', 'MarkerFaceColor', 'b')
plot(Y, Y, 'k-')
ylim([-0.01 0.8])
xlim([-0.01 2.6])
title('Linear Sum vs Multisensory Response')
xlabel('linear sum mean')
ylabel('multisensory response mean')
legend('Vis', 'Mech', 'Y=X', 'location', 'best')
hold off

suptitle(sprintf('Summary of all exps: Principle of Inverse Effectiveness raw values, no outliers'))