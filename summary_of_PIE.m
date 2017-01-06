%% Make a summary plot with subplots for each experiment. 

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

figure;

% subplot 1: area - lgst uni vs MSenh all
% get data to plot, by type of max unisensory
[unisensory, uniloc] = max(tadpole.area_avg(2:3,:),[],1); %get val and loc of max unisensory
area_V = [unisensory(uniloc==1); tadpole.area_avg(uniloc==1)]; % where max unisensory = V (1)
area_M = [unisensory(uniloc==2); tadpole.area_avg(uniloc==2)]; % where max unisensory = M (2)

% get values for plotting Y=X line
uni_mm = [floor(min(min([area_V(1,:) area_M(1,:)]))), ceil(max(max([area_V(1,:), area_M(1,:)])))];
multi_mm = [floor(min(min([area_V(2,:), area_M(2,:)]))), ceil(max(max([area_V(2,:), area_M(2,:)])))];
Y = (min([uni_mm, multi_mm])-1):1:(max(([uni_mm, multi_mm]))+1);

subplot(2,5,1)
hold on
plot(area_V(1,:), area_V(2,:), 'g*')
plot(area_M(1,:), area_M(2,:), 'm*')
plot(Y, Y, 'b-')
title('Area: largest uni vs MS enh')
xlabel('largest unisensory mean')
ylabel('MS enhancement mean')
legend('Vis', 'Mech', 'Y=X', 'location', 'best')
hold off

% subplot 2: area - lgst uni vs MSenh positive only
%%%%%%%%%%%%%%%%% SOMETHING IS WRONG HERE.
% get data to plot, by type of max unisensory, positive vals only
area_Vpos = [area_V(:, area_V(1,:) > 0)];
area_Mpos = [area_M(:, area_M(1,:) > 0)];
uni_mm_pos = [floor(min(min([area_Vpos(1,:) area_Mpos(1,:)]))), ceil(max(max([area_Vpos(1,:), area_Mpos(1,:)])))];
multi_mm_pos = [floor(min(min([area_Vpos(2,:), area_Mpos(2,:)]))), ceil(max(max([area_Vpos(2,:), area_Mpos(2,:)])))];
Y = (min([uni_mm_pos, multi_mm_pos])-1):1:(max(([uni_mm_pos, multi_mm_pos]))+1);

subplot(2,5,2)
hold on
plot(area_Vpos(1,:), area_Vpos(2,:), 'g*')
plot(area_Mpos(1,:), area_Mpos(2,:), 'm*')
plot(Y, Y, 'b-')
title('Area: largest uni vs MS enh, positive only')
xlabel('largest unisensory mean')
ylabel('MS enhancement mean')
legend('Vis', 'Mech', 'Y=X', 'location', 'best')
hold off

% subplot 3: area - lgst uni vs MSenh negative only 
%%%%%%%%%%%%%%%%% SOMETHING IS WRONG HERE.
% get data to plot, by type of max unisensory, negative vals only
area_Vneg = [area_V(:, area_V(1,:) < 0)];
area_Mneg = [area_M(:, area_M(1,:) < 0)];
uni_mm_neg  = [floor(min(min([area_Vneg(1,:) area_Mneg(1,:)]))), ceil(max(max([area_Vneg(1,:), area_Mneg(1,:)])))];
multi_mm_neg  = [floor(min(min([area_Vneg(2,:), area_Mneg(2,:)]))), ceil(max(max([area_Vneg(2,:), area_Mneg(2,:)])))];
Y = (min([uni_mm_neg , multi_mm_neg ])-1):1:(max(([uni_mm_neg , multi_mm_neg ]))+1);

subplot(2,5,3)
hold on
plot(area_Vneg(1,:), area_Vneg(2,:), 'g*')
plot(area_Mneg(1,:), area_Mneg(2,:), 'm*')
plot(Y, Y, 'b-')
title('Area: largest uni vs MS enh, negative only')
xlabel('largest unisensory mean')
ylabel('MS enhancement mean')
legend('Vis', 'Mech', 'Y=X', 'location', 'best')
hold off

% subplot 4: area - linear sum vs MS response
% get data to plot 
area_Vls = tadpole.area_avg(2,uniloc==1) + tadpole.area_avg(3,uniloc == 1);
area_Vm = tadpole.area_avg(1,uniloc==1);
area_Mls = tadpole.area_avg(2,uniloc==2) + tadpole.area_avg(3,uniloc == 2);
area_Mm = tadpole.area_avg(1,uniloc==2);

% get values for plotting Y=X line
uni_mm_ls = [floor(min(min([area_Vls(1,:) area_Mls(1,:)]))), ceil(max(max([area_Vls(1,:), area_Mls(1,:)])))];
multi_mm_m = [floor(min(min([area_Vm(1,:), area_Mm(1,:)]))), ceil(max(max([area_Vm(1,:), area_Mm(1,:)])))];
Y = (min([uni_mm_ls, multi_mm_m])-1):1:(max(([uni_mm_ls, multi_mm_m]))+1);

subplot(2,5,4)
hold on
plot(area_Vls(1,:), area_Vm(1,:), 'g*')
plot(area_Mls(1,:), area_Mm(1,:), 'm*')
plot(Y, Y, 'b-')
title('Area: linear sum vs multisensory')
xlabel('linear sum mean')
ylabel('multisensory mean')
legend('Vis', 'Mech', 'Y=X', 'location', 'best')
hold off

% subplot 5: area histogram of no stim trials
% get data for mean and stdev
stdev = std(cell2mat(tadpole.stim_vals_area{1,4}));
mean_val = mean(cell2mat(tadpole.stim_vals_area{1,4}));

hax = subplot(2,5,5);
hold on
%hax=axes;
hist(cell2mat(tadpole.stim_vals_area{1,4}),50)
line([(stdev+mean_val) (stdev+mean_val)], get(hax,'YLim'),'Color',[1 0 0])
line([(-stdev+mean_val) (-stdev+mean_val)], get(hax,'YLim'),'Color',[1 0 0])
line([(2*stdev+mean_val) (2*stdev+mean_val)], get(hax,'YLim'),'Color',[0 0 1])
line([(-2*stdev+mean_val) (-2*stdev+mean_val)], get(hax,'YLim'),'Color',[0 0 1])
line([mean_val mean_val], get(hax,'YLim'),'Color',[0 1 0])
title('All area measurements for no stimulus trials')
xlabel('area')
ylabel('counts')
hold off


% now by peak:
% subplot 6: peak - lgst uni vs MSenh all
% get data to plot, by type of max unisensory
[unisensory, uniloc] = max(tadpole.peak_avg(2:3,:),[],1); %get val and loc of max unisensory
peak_V = [unisensory(uniloc==1); tadpole.peak_avg(uniloc==1)]; % where max unisensory = V (1)
peak_M = [unisensory(uniloc==2); tadpole.peak_avg(uniloc==2)]; % where max unisensory = M (2)

% get values for plotting Y=X line
uni_mm = [min(min([peak_V(1,:) peak_M(1,:)])), max(max([peak_V(1,:), peak_M(1,:)]))];
multi_mm = [min(min([peak_V(2,:), peak_M(2,:)])), max(max([peak_V(2,:), peak_M(2,:)]))];
Y = (min([uni_mm, multi_mm])):0.01:(max(([uni_mm, multi_mm])));

subplot(2,5,6)
hold on
plot(peak_V(1,:), peak_V(2,:), 'g*')
plot(peak_M(1,:), peak_M(2,:), 'm*')
plot(Y, Y, 'b-')
title('peak: largest uni vs MS enh')
xlabel('largest unisensory mean')
ylabel('MS enhancement mean')
legend('Vis', 'Mech', 'Y=X', 'location', 'best')
hold off

% subplot 7: area - lgst uni vs MSenh positive only
%%%%%%%%%%%%%%%%% SOMETHING IS WRONG HERE.
% get data to plot, by type of max unisensory, positive vals only
peak_Vpos = [peak_V(:, peak_V(1,:) > 0)];
peak_Mpos = [peak_M(:, peak_M(1,:) > 0)];
uni_mm_pos = [min(min([peak_Vpos(1,:) peak_Mpos(1,:)])), max(max([peak_Vpos(1,:), peak_Mpos(1,:)]))];
multi_mm_pos = [min(min([peak_Vpos(2,:), peak_Mpos(2,:)])), max(max([peak_Vpos(2,:), peak_Mpos(2,:)]))];
Y = (min([uni_mm_pos, multi_mm_pos])):0.01:(max(([uni_mm_pos, multi_mm_pos])));

subplot(2,5,7)
hold on
plot(peak_Vpos(1,:), peak_Vpos(2,:), 'g*')
plot(peak_Mpos(1,:), peak_Mpos(2,:), 'm*')
plot(Y, Y, 'b-')
title('peak: largest uni vs MS enh, positive only')
xlabel('largest unisensory mean')
ylabel('MS enhancement mean')
legend('Vis', 'Mech', 'Y=X', 'location', 'best')
hold off

% subplot 8: area - lgst uni vs MSenh negative only 
%%%%%%%%%%%%%%%%% SOMETHING IS WRONG HERE.
% get data to plot, by type of max unisensory, negative vals only
peak_Vneg = [peak_V(:, peak_V(1,:) < 0)];
peak_Mneg = [peak_M(:, peak_M(1,:) < 0)];
uni_mm_neg  = [min(min([peak_Vneg(1,:) peak_Mneg(1,:)])), max(max([peak_Vneg(1,:), peak_Mneg(1,:)]))];
multi_mm_neg  = [min(min([peak_Vneg(2,:), peak_Mneg(2,:)])), max(max([peak_Vneg(2,:), peak_Mneg(2,:)]))];
Y = (min([uni_mm_neg , multi_mm_neg ])):0.01:(max(([uni_mm_neg , multi_mm_neg ])));

subplot(2,5,8)
hold on
plot(peak_Vneg(1,:), peak_Vneg(2,:), 'g*')
plot(peak_Mneg(1,:), peak_Mneg(2,:), 'm*')
plot(Y, Y, 'b-')
title('peak: largest uni vs MS enh, negative only')
xlabel('largest unisensory mean')
ylabel('MS enhancement mean')
legend('Vis', 'Mech', 'Y=X', 'location', 'best')
hold off

% subplot 9: area - linear sum vs MS response
% get data to plot 
peak_Vls = tadpole.peak_avg(2,uniloc==1) + tadpole.peak_avg(3,uniloc == 1);
peak_Vm = tadpole.peak_avg(1,uniloc==1);
peak_Mls = tadpole.peak_avg(2,uniloc==2) + tadpole.peak_avg(3,uniloc == 2);
peak_Mm = tadpole.peak_avg(1,uniloc==2);

% get values for plotting Y=X line
uni_mm_ls = [min(min([peak_Vls(1,:) peak_Mls(1,:)])), max(max([peak_Vls(1,:), peak_Mls(1,:)]))];
multi_mm_m = [min(min([peak_Vm(1,:), peak_Mm(1,:)])), max(max([peak_Vm(1,:), peak_Mm(1,:)]))];
Y = (min([uni_mm_ls, multi_mm_m])):0.01:(max(([uni_mm_ls, multi_mm_m])));

subplot(2,5,9)
hold on
plot(peak_Vls(1,:), peak_Vm(1,:), 'g*')
plot(peak_Mls(1,:), peak_Mm(1,:), 'm*')
plot(Y, Y, 'b-')
title('peak: linear sum vs multisensory')
xlabel('linear sum mean')
ylabel('multisensory mean')
legend('Vis', 'Mech', 'Y=X', 'location', 'best')
hold off

% subplot 5: area histogram of no stim trials
% get data for mean and stdev
stdev = std(cell2mat(tadpole.stim_vals_meanpeak{1,4}));
mean_val = mean(cell2mat(tadpole.stim_vals_meanpeak{1,4}));

hax = subplot(2,5,10);
hold on

hist(cell2mat(tadpole.stim_vals_meanpeak{1,4}),50)
line([(stdev+mean_val) (stdev+mean_val)], get(hax,'YLim'),'Color',[1 0 0])
line([(-stdev+mean_val) (-stdev+mean_val)], get(hax,'YLim'),'Color',[1 0 0])
line([(2*stdev+mean_val) (2*stdev+mean_val)], get(hax,'YLim'),'Color',[0 0 1])
line([(-2*stdev+mean_val) (-2*stdev+mean_val)], get(hax,'YLim'),'Color',[0 0 1])
line([mean_val mean_val], get(hax,'YLim'),'Color',[0 1 0])
title('All peak measurements for no stimulus trials')
xlabel('peak')
ylabel('counts')
hold off

%main title
suptitle(sprintf([matFiles(k).name, ' Principle of Inverse Effectiveness raw values']));
%suptitle(sprintf('Exp 5 blocks 6-9 Principle of Inverse Effectiveness raw values'))
%fig_filename=(sprintf([myFolder, 'figures/', matFiles(k).name, 'PIE_summary_rawvals']));
%fig_filename=(sprintf([myFolder, 'figures/', 'exp20' 'PIE_summary_rawvals']));
%saveas(gcf,fig_filename,'png');
clear('tadpole');
end