%% Final PIE Paper Calcium traces (with adjustments for bleaching in all displays)

%% Exp 5 ROI 51 (MSenh)
% exp 5 is in tadpole{1,7}

% plot all traces
%figure;
t = 7;
% identify trials that equal each stim type
clear('stimorder_uniques', 'stimlist')
stimorder_uniques = unique(tadpole{1,t}.stimorder)
for s = 1:length(stimorder_uniques)
    stimlist(s,:) = find(tadpole{1,t}.stimorder == stimorder_uniques(s));
end
color_list = [1 1 0; 1 0 1; 0 1 1; 1 0 0; 0 1 0; 0 0 1; 0 0 0; 0.5 0 0.5; 0.5 0.5 0.5; 0.5 0.5 0; 0 0.5 0.5; 0 0 0.5]     
r = 51;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% averages for each modality
% generate matrix of traces I need
clear('filteredMS', 'filteredM', 'filteredV', 'avgtraceMS', 'avgtraceM', 'avgtraceV')
% multisensory
for i = 1:size(stimlist,2)
    filteredMS(i,:) = tadpole{1,t}.filtered{r, stimlist(1,i)};
end
avgtraceMS = mean(filteredMS);
% visual
for i = 1:size(stimlist,2)
    filteredV(i,:) = tadpole{1,t}.filtered{r, stimlist(2,i)};
end
avgtraceV = mean(filteredV);
%mechanosensory
for i = 1:size(stimlist,2)
    filteredM(i,:) = tadpole{1,t}.filtered{r, stimlist(3,i)};
end
avgtraceM = mean(filteredM);
% no stim
for i = 1:size(stimlist,2)
    filteredNS(i,:) = tadpole{1,t}.filtered{r, stimlist(4,i)};
end
avgtraceNS = mean(filteredNS);

% figure out slope of decay from dye bleaching using no stim trials
% fit a single polynomial (e.g. a line) to each filtered dataset
X_vals = 1:1:160;
for i = 1:size(filteredNS, 1) % over each presentation
    [p(i,:), S(i), mu(i,:)] = polyfit(X_vals, filteredNS(i, :), 1)
    %p(i) = polyfit(X_vals, filteredNS(i, :), 1)
end

avg_slope = mean(p([2:5,7,8],1))
stddev_slope = std(p(:,1))
add_on = -avg_slope * X_vals *0.03

% correct slope values by add_on - avg traces
avgtraceM_adj = avgtraceM + add_on
avgtraceV_adj = avgtraceV + add_on
avgtraceMS_adj = avgtraceMS + add_on
linsum_adj = avgtraceM_adj + avgtraceV_adj
% correct slope vals by add_on - individual traces
for i = 1:size(filteredV, 1)
    filteredV_adj(i,:) = filteredV(i,:) + add_on;
end
for i = 1:size(filteredV, 1)
    filteredM_adj(i,:) = filteredM(i,:) + add_on;
end
for i = 1:size(filteredV, 1)
    filteredMS_adj(i,:) = filteredMS(i,:) + add_on;
end

% plot avg traces for Ms, M, V
figure;
hold on
plot(avgtraceMS_adj, 'color', [0.5 0 0.5])
plot(avgtraceV_adj, 'color', 'b')
plot(avgtraceM_adj, 'color', 'r')
plot(linsum_adj, 'color', 'k')
hold off
title(sprintf('mean responses exp %d roi %d bleach adj', tadpole{1,t}.expnum, r))
xlabel('time (s)')
ylabel('\DeltaF/F_{0}')
ax=gca;
xsize = length(tadpole{1,t}.filtered{1,1});
y_min = min([min(avgtraceM_adj), min(avgtraceV_adj), min(avgtraceMS_adj), min(linsum_adj)])
y_max = max([max(avgtraceM_adj), max(avgtraceV_adj), max(avgtraceMS_adj), max(linsum_adj)])
ax.XTick = [0, xsize/7, (xsize/7)*2, (xsize/7)*3, (xsize/7)*4, (xsize/7)*5, (xsize/7)*6, (xsize/7)*7];
ax.XTickLabel = {'0','1', '2', '3', '4', '5', '6', '7'};
axis([2 xsize-1 y_min y_max])
% print to file (color!!)
clear('fig_filename')
fig_filename=sprintf('F:/Calcium_Imaging_Analysis/tadpoles_byexp/figures/smoothedROIs/paper_examples/adjusted/exp%droi%dsmoothed_averageall.eps', tadpole{1,t}.expnum, r)
% save as eps for illustrator (in color)
saveas(gcf,fig_filename,'psc2');
% save as PNG for easy perusal
clear('fig_filename')
fig_filename=sprintf('F:/Calcium_Imaging_Analysis/tadpoles_byexp/figures/smoothedROIs/paper_examples/adjusted/exp%droi%dsmoothed_averageall.png', tadpole{1,t}.expnum, r)
saveas(gcf,fig_filename,'png');

% exact plots desired for paper:
tracesMS = [1, 12, 4] %(Y, C, R)
tracesV = [4, 6, 7] %(R, blue, black)
tracesM = [1, 5, 7] %(Y, gr, black)

% axis dims for individual plots
y_minMS = min([min(filteredMS_adj(tracesMS(1,1),:)), min(filteredMS_adj(tracesMS(1,2),:)), min(filteredMS_adj(tracesMS(1,3),:))])
y_maxMS = max([max(filteredMS_adj(tracesMS(1,1),:)), max(filteredMS_adj(tracesMS(1,2),:)), max(filteredMS_adj(tracesMS(1,3),:))])
y_minM = min([min(filteredM_adj(tracesM(1,1),:)), min(filteredM_adj(tracesM(1,2),:)), min(filteredM_adj(tracesM(1,3),:))])
y_maxM = max([max(filteredM_adj(tracesM(1,1),:)), max(filteredM_adj(tracesM(1,2),:)), max(filteredM_adj(tracesM(1,3),:))])
y_minV = min([min(filteredV_adj(tracesV(1,1),:)), min(filteredV_adj(tracesV(1,2),:)), min(filteredV_adj(tracesV(1,3),:))])
y_maxV = max([max(filteredV_adj(tracesV(1,1),:)), max(filteredV_adj(tracesV(1,2),:)), max(filteredV_adj(tracesV(1,3),:))])
y_min = min([y_minMS, y_minM, y_minV])
y_max = max([y_maxMS, y_maxM, y_maxV])
xsize = size(filteredMS_adj,2);

% multisensory (purple)
figure;
hold on
for k = 1:length(tracesMS)
    plot(filteredMS_adj(tracesMS(k),:), 'color', [0.5 0 0.5])
end
title(sprintf('Exp %d ROI %d smoothed multisensory', tadpole{1,t}.expnum, r))
xlabel('time (s)')
ylabel('\DeltaF/F_{0}')
ax=gca;

y_min = min([min(filteredMS_adj(tracesMS(1,1),:)), min(filteredMS_adj(tracesMS(1,2),:)), min(filteredMS_adj(tracesMS(1,3),:))])
y_max = max([max(filteredMS_adj(tracesMS(1,1),:)), max(filteredMS_adj(tracesMS(1,2),:)), max(filteredMS_adj(tracesMS(1,3),:))])
ax.XTick = [0, xsize/7, (xsize/7)*2, (xsize/7)*3, (xsize/7)*4, (xsize/7)*5, (xsize/7)*6, (xsize/7)*7];
ax.XTickLabel = {'0','1', '2', '3', '4', '5', '6', '7'};
hold off
axis([2 xsize y_min y_max])

fig_filename=sprintf('F:/Calcium_Imaging_Analysis/tadpoles_byexp/figures/smoothedROIs/paper_examples/adjusted/exp%droi%dsmoothed_multiselected.eps', tadpole{1,t}.expnum, r)
saveas(gcf,fig_filename,'psc2');

% visual (blue)
figure;
hold on
for k = 1:length(tracesV)
    plot(filteredV_adj(tracesV(k),:), 'color', 'b')
end
title(sprintf('Exp %d ROI %d smoothed visual', tadpole{1,t}.expnum, r))
xlabel('time (s)')
ylabel('\DeltaF/F_{0}')
ax=gca;
xsize = size(filteredMS_adj,2);
y_min = min([min(filteredMS_adj(tracesMS(1,1),:)), min(filteredMS_adj(tracesMS(1,2),:)), min(filteredMS_adj(tracesMS(1,3),:))])
y_max = max([max(filteredMS_adj(tracesMS(1,1),:)), max(filteredMS_adj(tracesMS(1,2),:)), max(filteredMS_adj(tracesMS(1,3),:))])
ax.XTick = [0, xsize/7, (xsize/7)*2, (xsize/7)*3, (xsize/7)*4, (xsize/7)*5, (xsize/7)*6, (xsize/7)*7];
ax.XTickLabel = {'0','1', '2', '3', '4', '5', '6', '7'};
hold off
axis([2 xsize y_min y_max])
fig_filename=sprintf('F:/Calcium_Imaging_Analysis/tadpoles_byexp/figures/smoothedROIs/paper_examples/adjusted/exp%droi%dsmoothed_visualselected.eps', tadpole{1,t}.expnum, r)
saveas(gcf,fig_filename,'psc2');

% mechanosensory (red)
figure;
hold on
for k = 1:length(tracesM)
    plot(filteredM_adj(tracesV(k),:), 'color', 'r')
end
title(sprintf('Exp %d ROI %d smoothed mech', tadpole{1,t}.expnum, r))
xlabel('time (s)')
ylabel('\DeltaF/F_{0}')
ax=gca;
ax.XTick = [0, xsize/7, (xsize/7)*2, (xsize/7)*3, (xsize/7)*4, (xsize/7)*5, (xsize/7)*6, (xsize/7)*7];
ax.XTickLabel = {'0','1', '2', '3', '4', '5', '6', '7'};
hold off
axis([2 xsize y_min y_max])
fig_filename=sprintf('F:/Calcium_Imaging_Analysis/tadpoles_byexp/figures/smoothedROIs/paper_examples/adjusted/exp%droi%dsmoothed_mechselected.eps', tadpole{1,t}.expnum, r)
saveas(gcf,fig_filename,'psc2');

%% Exp 5 ROI 151 (MSenh)
% plot all traces
t = 7;
% identify trials that equal each stim type
clear('stimorder_uniques', 'stimlist')
stimorder_uniques = unique(tadpole{1,t}.stimorder)
for s = 1:length(stimorder_uniques)
    stimlist(s,:) = find(tadpole{1,t}.stimorder == stimorder_uniques(s));
end
color_list = [1 1 0; 1 0 1; 0 1 1; 1 0 0; 0 1 0; 0 0 1; 0 0 0; 0.5 0 0.5; 0.5 0.5 0.5; 0.5 0.5 0; 0 0.5 0.5; 0 0 0.5]     
r = 151;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% averages for each modality
% generate matrix of traces I need
clear('filteredMS', 'filteredM', 'filteredV', 'avgtraceMS', 'avgtraceM', 'avgtraceV', 'linsum_adj', ...
    'filteredMS_adj', 'filteredM_adj', 'filteredV_adj', 'avgtraceMS_adj', 'avgtraceM_adj', 'avgtraceV_adj')
% multisensory
for i = 1:size(stimlist,2)
    filteredMS(i,:) = tadpole{1,t}.filtered{r, stimlist(1,i)};
end
avgtraceMS = mean(filteredMS);
% visual
for i = 1:size(stimlist,2)
    filteredV(i,:) = tadpole{1,t}.filtered{r, stimlist(2,i)};
end
avgtraceV = mean(filteredV);
%mechanosensory
for i = 1:size(stimlist,2)
    filteredM(i,:) = tadpole{1,t}.filtered{r, stimlist(3,i)};
end
avgtraceM = mean(filteredM);
% no stim
for i = 1:size(stimlist,2)
    filteredNS(i,:) = tadpole{1,t}.filtered{r, stimlist(4,i)};
end
avgtraceNS = mean(filteredNS);

% figure out slope of decay from dye bleaching using no stim trials
% fit a single polynomial (e.g. a line) to each filtered dataset
X_vals = 1:1:160;
for i = 1:size(filteredNS, 1) % over each presentation
    [p(i,:), S(i), mu(i,:)] = polyfit(X_vals, filteredNS(i, :), 1)
    %p(i) = polyfit(X_vals, filteredNS(i, :), 1)
end

avg_slope = mean(p([2:5,7,8],1))
stddev_slope = std(p(:,1))
add_on = -avg_slope * X_vals *0.15

% correct slope values by add_on - avg traces
avgtraceM_adj = avgtraceM + add_on
avgtraceV_adj = avgtraceV + add_on
avgtraceMS_adj = avgtraceMS + add_on
linsum_adj = avgtraceM_adj + avgtraceV_adj
% correct slope vals by add_on - individual traces
for i = 1:size(filteredV, 1)
    filteredV_adj(i,:) = filteredV(i,:) + add_on;
end
for i = 1:size(filteredV, 1)
    filteredM_adj(i,:) = filteredM(i,:) + add_on;
end
for i = 1:size(filteredV, 1)
    filteredMS_adj(i,:) = filteredMS(i,:) + add_on;
end

% plot avg traces for Ms, M, V
figure;
hold on
plot(avgtraceMS_adj, 'color', [0.5 0 0.5])
plot(avgtraceV_adj, 'color', 'b')
plot(avgtraceM_adj, 'color', 'r')
plot(linsum_adj, 'color', 'k')
hold off
title(sprintf('mean responses exp %d roi %d bleach adj', tadpole{1,t}.expnum, r))
xlabel('time (s)')
ylabel('\DeltaF/F_{0}')
ax=gca;
xsize = length(tadpole{1,t}.filtered{1,1});
y_min = min([min(avgtraceM_adj), min(avgtraceV_adj), min(avgtraceMS_adj), min(linsum_adj)])
y_max = max([max(avgtraceM_adj), max(avgtraceV_adj), max(avgtraceMS_adj), max(linsum_adj)])
ax.XTick = [0, xsize/7, (xsize/7)*2, (xsize/7)*3, (xsize/7)*4, (xsize/7)*5, (xsize/7)*6, (xsize/7)*7];
ax.XTickLabel = {'0','1', '2', '3', '4', '5', '6', '7'};
axis([2 xsize-1 y_min y_max])
% print to file (color!!)
clear('fig_filename')
fig_filename=sprintf('F:/Calcium_Imaging_Analysis/tadpoles_byexp/figures/smoothedROIs/paper_examples/adjusted/exp%droi%dsmoothed_averageall.eps', tadpole{1,t}.expnum, r)
% save as eps for illustrator (in color)
saveas(gcf,fig_filename,'psc2');
% save as PNG for easy perusal
clear('fig_filename')
fig_filename=sprintf('F:/Calcium_Imaging_Analysis/tadpoles_byexp/figures/smoothedROIs/paper_examples/adjusted/exp%droi%dsmoothed_averageall.png', tadpole{1,t}.expnum, r)
saveas(gcf,fig_filename,'png');

% exact plots desired for paper:
tracesMS = [1, 4, 6] %(Y, C, R)
tracesV = [2, 3, 5] %(R, blue, black)
tracesM = [1, 2, 4] %(Y, gr, black)

% axis dims for individual plots
y_minMS = min([min(filteredMS_adj(tracesMS(1,1),:)), min(filteredMS_adj(tracesMS(1,2),:)), min(filteredMS_adj(tracesMS(1,3),:))])
y_maxMS = max([max(filteredMS_adj(tracesMS(1,1),:)), max(filteredMS_adj(tracesMS(1,2),:)), max(filteredMS_adj(tracesMS(1,3),:))])
y_minM = min([min(filteredM_adj(tracesM(1,1),:)), min(filteredM_adj(tracesM(1,2),:)), min(filteredM_adj(tracesM(1,3),:))])
y_maxM = max([max(filteredM_adj(tracesM(1,1),:)), max(filteredM_adj(tracesM(1,2),:)), max(filteredM_adj(tracesM(1,3),:))])
y_minV = min([min(filteredV_adj(tracesV(1,1),:)), min(filteredV_adj(tracesV(1,2),:)), min(filteredV_adj(tracesV(1,3),:))])
y_maxV = max([max(filteredV_adj(tracesV(1,1),:)), max(filteredV_adj(tracesV(1,2),:)), max(filteredV_adj(tracesV(1,3),:))])
y_min = min([y_minMS, y_minM, y_minV])
y_max = max([y_maxMS, y_maxM, y_maxV])
xsize = size(filteredMS_adj,2);

% multisensory (purple)
figure;
hold on
for k = 1:length(tracesMS)
    plot(filteredMS_adj(tracesMS(k),:), 'color', [0.5 0 0.5])
end
title(sprintf('Exp %d ROI %d smoothed multisensory', tadpole{1,t}.expnum, r))
xlabel('time (s)')
ylabel('\DeltaF/F_{0}')
ax=gca;

y_min = min([min(filteredMS_adj(tracesMS(1,1),:)), min(filteredMS_adj(tracesMS(1,2),:)), min(filteredMS_adj(tracesMS(1,3),:))])
y_max = max([max(filteredMS_adj(tracesMS(1,1),:)), max(filteredMS_adj(tracesMS(1,2),:)), max(filteredMS_adj(tracesMS(1,3),:))])
ax.XTick = [0, xsize/7, (xsize/7)*2, (xsize/7)*3, (xsize/7)*4, (xsize/7)*5, (xsize/7)*6, (xsize/7)*7];
ax.XTickLabel = {'0','1', '2', '3', '4', '5', '6', '7'};
hold off
axis([2 xsize y_min y_max])

fig_filename=sprintf('F:/Calcium_Imaging_Analysis/tadpoles_byexp/figures/smoothedROIs/paper_examples/adjusted/exp%droi%dsmoothed_multiselected.eps', tadpole{1,t}.expnum, r)
saveas(gcf,fig_filename,'psc2');

% visual (blue)
figure;
hold on
for k = 1:length(tracesV)
    plot(filteredV_adj(tracesV(k),:), 'color', 'b')
end
title(sprintf('Exp %d ROI %d smoothed visual', tadpole{1,t}.expnum, r))
xlabel('time (s)')
ylabel('\DeltaF/F_{0}')
ax=gca;
xsize = size(filteredMS_adj,2);
y_min = min([min(filteredMS_adj(tracesMS(1,1),:)), min(filteredMS_adj(tracesMS(1,2),:)), min(filteredMS_adj(tracesMS(1,3),:))])
y_max = max([max(filteredMS_adj(tracesMS(1,1),:)), max(filteredMS_adj(tracesMS(1,2),:)), max(filteredMS_adj(tracesMS(1,3),:))])
ax.XTick = [0, xsize/7, (xsize/7)*2, (xsize/7)*3, (xsize/7)*4, (xsize/7)*5, (xsize/7)*6, (xsize/7)*7];
ax.XTickLabel = {'0','1', '2', '3', '4', '5', '6', '7'};
hold off
axis([2 xsize y_min y_max])
fig_filename=sprintf('F:/Calcium_Imaging_Analysis/tadpoles_byexp/figures/smoothedROIs/paper_examples/adjusted/exp%droi%dsmoothed_visualselected.eps', tadpole{1,t}.expnum, r)
saveas(gcf,fig_filename,'psc2');

% mechanosensory (red)
figure;
hold on
for k = 1:length(tracesM)
    plot(filteredM_adj(tracesV(k),:), 'color', 'r')
end
title(sprintf('Exp %d ROI %d smoothed mech', tadpole{1,t}.expnum, r))
xlabel('time (s)')
ylabel('\DeltaF/F_{0}')
ax=gca;
ax.XTick = [0, xsize/7, (xsize/7)*2, (xsize/7)*3, (xsize/7)*4, (xsize/7)*5, (xsize/7)*6, (xsize/7)*7];
ax.XTickLabel = {'0','1', '2', '3', '4', '5', '6', '7'};
hold off
axis([2 xsize y_min y_max])
fig_filename=sprintf('F:/Calcium_Imaging_Analysis/tadpoles_byexp/figures/smoothedROIs/paper_examples/adjusted/exp%droi%dsmoothed_mechselected.eps', tadpole{1,t}.expnum, r)
saveas(gcf,fig_filename,'psc2');

%% Exp 6 ROI 45
t = 8;
% identify trials that equal each stim type
clear('stimorder_uniques', 'stimlist')
stimorder_uniques = unique(tadpole{1,t}.stimorder)
for s = 1:length(stimorder_uniques)
    stimlist(s,:) = find(tadpole{1,t}.stimorder == stimorder_uniques(s));
end
color_list = [1 1 0; 1 0 1; 0 1 1; 1 0 0; 0 1 0; 0 0 1; 0 0 0; 0.5 0 0.5; 0.5 0.5 0.5; 0.5 0.5 0; 0 0.5 0.5; 0 0 0.5]     
r = 45;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% averages for each modality
% generate matrix of traces I need
clear('filteredMS', 'filteredM', 'filteredV', 'filteredNS', 'avgtraceMS', 'avgtraceM', 'avgtraceV', 'linsum_adj', ...
    'filteredMS_adj', 'filteredM_adj', 'filteredV_adj', 'avgtraceMS_adj', 'avgtraceM_adj', 'avgtraceV_adj')
% multisensory
for i = 1:6
    filteredMS(i,:) = tadpole{1,t}.filtered{r, stimlist(1,i)};
end
avgtraceMS = mean(filteredMS);
% visual
for i = 1:6
    filteredV(i,:) = tadpole{1,t}.filtered{r, stimlist(2,i)};
end
avgtraceV = mean(filteredV);
%mechanosensory
for i = 1:6
    filteredM(i,:) = tadpole{1,t}.filtered{r, stimlist(3,i)};
end
avgtraceM = mean(filteredM);
% no stim
for i = 1:6
    filteredNS(i,:) = tadpole{1,t}.filtered{r, stimlist(4,i)};
end
avgtraceNS = mean(filteredNS);

% figure out slope of decay from dye bleaching using no stim trials
% fit a single polynomial (e.g. a line) to each filtered dataset
X_vals = 1:1:220;
for i = 1:size(filteredNS, 1) % over each presentation
    [p(i,:), S(i), mu(i,:)] = polyfit(X_vals, filteredNS(i, :), 1)
    %p(i) = polyfit(X_vals, filteredNS(i, :), 1)
end

avg_slope = mean(p([2:5,7,8],1))
stddev_slope = std(p(:,1))
add_on = -avg_slope * X_vals *0.005

% correct slope values by add_on - avg traces
avgtraceM_adj = avgtraceM + add_on
avgtraceV_adj = avgtraceV + add_on
avgtraceMS_adj = avgtraceMS + add_on
linsum_adj = avgtraceM_adj + avgtraceV_adj
% correct slope vals by add_on - individual traces
for i = 1:size(filteredV, 1)
    filteredV_adj(i,:) = filteredV(i,:) + add_on;
end
for i = 1:size(filteredV, 1)
    filteredM_adj(i,:) = filteredM(i,:) + add_on;
end
for i = 1:size(filteredV, 1)
    filteredMS_adj(i,:) = filteredMS(i,:) + add_on;
end

% plot avg traces for Ms, M, V
figure;
hold on
plot(avgtraceMS_adj, 'color', [0.5 0 0.5])
plot(avgtraceV_adj, 'color', 'b')
plot(avgtraceM_adj, 'color', 'r')
plot(linsum_adj, 'color', 'k')
hold off
title(sprintf('mean responses exp %d roi %d bleach adj', tadpole{1,t}.expnum, r))
xlabel('time (s)')
ylabel('\DeltaF/F_{0}')
ax=gca;
xsize = length(tadpole{1,t}.filtered{1,1});
y_min = min([min(avgtraceM_adj), min(avgtraceV_adj), min(avgtraceMS_adj), min(linsum_adj)])
y_max = max([max(avgtraceM_adj), max(avgtraceV_adj), max(avgtraceMS_adj), max(linsum_adj)])
ax.XTick = [0, xsize/7, (xsize/7)*2, (xsize/7)*3, (xsize/7)*4, (xsize/7)*5, (xsize/7)*6, (xsize/7)*7];
ax.XTickLabel = {'0','1', '2', '3', '4', '5', '6', '7'};
axis([2 xsize-1 y_min y_max])
% print to file (color!!)
clear('fig_filename')
fig_filename=sprintf('F:/Calcium_Imaging_Analysis/tadpoles_byexp/figures/smoothedROIs/paper_examples/adjusted/exp%droi%dsmoothed_averageall.eps', tadpole{1,t}.expnum, r)
% save as eps for illustrator (in color)
saveas(gcf,fig_filename,'psc2');
% save as PNG for easy perusal
clear('fig_filename')
fig_filename=sprintf('F:/Calcium_Imaging_Analysis/tadpoles_byexp/figures/smoothedROIs/paper_examples/adjusted/exp%droi%dsmoothed_averageall.png', tadpole{1,t}.expnum, r)
saveas(gcf,fig_filename,'png');

% exact plots desired for paper:
tracesMS = [2, 4, 6] %
tracesV = [2, 4, 5] %
tracesM = [1, 5, 6] %

% axis dims for individual plots
y_minMS = min([min(filteredMS_adj(tracesMS(1,1),:)), min(filteredMS_adj(tracesMS(1,2),:)), min(filteredMS_adj(tracesMS(1,3),:))])
y_maxMS = max([max(filteredMS_adj(tracesMS(1,1),:)), max(filteredMS_adj(tracesMS(1,2),:)), max(filteredMS_adj(tracesMS(1,3),:))])
y_minM = min([min(filteredM_adj(tracesM(1,1),:)), min(filteredM_adj(tracesM(1,2),:)), min(filteredM_adj(tracesM(1,3),:))])
y_maxM = max([max(filteredM_adj(tracesM(1,1),:)), max(filteredM_adj(tracesM(1,2),:)), max(filteredM_adj(tracesM(1,3),:))])
y_minV = min([min(filteredV_adj(tracesV(1,1),:)), min(filteredV_adj(tracesV(1,2),:)), min(filteredV_adj(tracesV(1,3),:))])
y_maxV = max([max(filteredV_adj(tracesV(1,1),:)), max(filteredV_adj(tracesV(1,2),:)), max(filteredV_adj(tracesV(1,3),:))])
y_min = min([y_minMS, y_minM, y_minV])
y_max = max([y_maxMS, y_maxM, y_maxV])
xsize = size(filteredMS_adj,2);

% multisensory (purple)
figure;
hold on
for k = 1:length(tracesMS)
    plot(filteredMS_adj(tracesMS(k),:), 'color', [0.5 0 0.5])
end
title(sprintf('Exp %d ROI %d smoothed multisensory', tadpole{1,t}.expnum, r))
xlabel('time (s)')
ylabel('\DeltaF/F_{0}')
ax=gca;
ax.XTick = [0, xsize/7, (xsize/7)*2, (xsize/7)*3, (xsize/7)*4, (xsize/7)*5, (xsize/7)*6, (xsize/7)*7];
ax.XTickLabel = {'0','1', '2', '3', '4', '5', '6', '7'};
hold off
axis([2 xsize-1 y_min y_max])

fig_filename=sprintf('F:/Calcium_Imaging_Analysis/tadpoles_byexp/figures/smoothedROIs/paper_examples/adjusted/exp%droi%dsmoothed_multiselected.eps', tadpole{1,t}.expnum, r)
saveas(gcf,fig_filename,'psc2');

% visual (blue)
figure;
hold on
for k = 1:length(tracesV)
    plot(filteredV_adj(tracesV(k),:), 'color', 'b')
end
title(sprintf('Exp %d ROI %d smoothed visual', tadpole{1,t}.expnum, r))
xlabel('time (s)')
ylabel('\DeltaF/F_{0}')
ax=gca;
ax.XTick = [0, xsize/7, (xsize/7)*2, (xsize/7)*3, (xsize/7)*4, (xsize/7)*5, (xsize/7)*6, (xsize/7)*7];
ax.XTickLabel = {'0','1', '2', '3', '4', '5', '6', '7'};
hold off
axis([2 xsize-1 y_min y_max])
fig_filename=sprintf('F:/Calcium_Imaging_Analysis/tadpoles_byexp/figures/smoothedROIs/paper_examples/adjusted/exp%droi%dsmoothed_visualselected.eps', tadpole{1,t}.expnum, r)
saveas(gcf,fig_filename,'psc2');

% mechanosensory (red)
figure;
hold on
for k = 1:length(tracesM)
    plot(filteredM_adj(tracesV(k),:), 'color', 'r')
end
title(sprintf('Exp %d ROI %d smoothed mech', tadpole{1,t}.expnum, r))
xlabel('time (s)')
ylabel('\DeltaF/F_{0}')
ax=gca;
ax.XTick = [0, xsize/7, (xsize/7)*2, (xsize/7)*3, (xsize/7)*4, (xsize/7)*5, (xsize/7)*6, (xsize/7)*7];
ax.XTickLabel = {'0','1', '2', '3', '4', '5', '6', '7'};
hold off
axis([2 xsize-1 y_min y_max])
fig_filename=sprintf('F:/Calcium_Imaging_Analysis/tadpoles_byexp/figures/smoothedROIs/paper_examples/adjusted/exp%droi%dsmoothed_mechselected.eps', tadpole{1,t}.expnum, r)
saveas(gcf,fig_filename,'psc2');

%% Exp 8 ROI 19
t = 10;
% identify trials that equal each stim type
clear('stimorder_uniques', 'stimlist')
stimorder_uniques = unique(tadpole{1,t}.stimorder)
for s = 1:length(stimorder_uniques)
    stimlist(s,:) = find(tadpole{1,t}.stimorder == stimorder_uniques(s));
end
color_list = [1 1 0; 1 0 1; 0 1 1; 1 0 0; 0 1 0; 0 0 1; 0 0 0; 0.5 0 0.5; 0.5 0.5 0.5; 0.5 0.5 0; 0 0.5 0.5; 0 0 0.5]     
r = 19;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% averages for each modality
% generate matrix of traces I need
clear('filteredMS', 'filteredM', 'filteredV', 'filteredNS', 'avgtraceMS', 'avgtraceM', 'avgtraceV', 'avgtraceNS', 'linsum_adj', ...
    'filteredMS_adj', 'filteredM_adj', 'filteredV_adj', 'avgtraceMS_adj', 'avgtraceM_adj', 'avgtraceV_adj')
% multisensory
for i = 4:size(stimlist,2)
    filteredMS(i,:) = tadpole{1,t}.filtered{r, stimlist(1,i)};
end
avgtraceMS = mean(filteredMS);
% visual
for i = 4:size(stimlist,2)
    filteredV(i,:) = tadpole{1,t}.filtered{r, stimlist(2,i)};
end
avgtraceV = mean(filteredV);
%mechanosensory
for i = 4:size(stimlist,2)
    filteredM(i,:) = tadpole{1,t}.filtered{r, stimlist(3,i)};
end
avgtraceM = mean(filteredM);
% no stim
for i = 4:size(stimlist,2)
    filteredNS(i,:) = tadpole{1,t}.filtered{r, stimlist(4,i)};
end
avgtraceNS = mean(filteredNS);

% figure out slope of decay from dye bleaching using no stim trials
% fit a single polynomial (e.g. a line) to each filtered dataset
X_vals = 1:1:160;
for i = 1:size(filteredNS, 1) % over each presentation
    [p(i,:), S(i), mu(i,:)] = polyfit(X_vals, filteredNS(i, :), 1)
    %p(i) = polyfit(X_vals, filteredNS(i, :), 1)
end

avg_slope = mean(p([2:5,7,8],1))
stddev_slope = std(p(:,1))
add_on = -avg_slope * X_vals *0.01

% correct slope values by add_on - avg traces
avgtraceM_adj = avgtraceM + add_on
avgtraceV_adj = avgtraceV + add_on
avgtraceMS_adj = avgtraceMS + add_on
linsum_adj = avgtraceM_adj + avgtraceV_adj
% correct slope vals by add_on - individual traces
for i = 1:size(filteredV, 1)
    filteredV_adj(i,:) = filteredV(i,:) + add_on;
end
for i = 1:size(filteredV, 1)
    filteredM_adj(i,:) = filteredM(i,:) + add_on;
end
for i = 1:size(filteredV, 1)
    filteredMS_adj(i,:) = filteredMS(i,:) + add_on;
end

% plot avg traces for Ms, M, V
figure;
hold on
plot(avgtraceMS_adj, 'color', [0.5 0 0.5])
plot(avgtraceV_adj, 'color', 'b')
plot(avgtraceM_adj, 'color', 'r')
plot(linsum_adj, 'color', 'k')
hold off
title(sprintf('mean responses exp %d roi %d bleach adj', tadpole{1,t}.expnum, r))
xlabel('time (s)')
ylabel('\DeltaF/F_{0}')
ax=gca;
xsize = length(tadpole{1,t}.filtered{1,1});
y_min = min([min(avgtraceM_adj), min(avgtraceV_adj), min(avgtraceMS_adj), min(linsum_adj)])
y_max = max([max(avgtraceM_adj), max(avgtraceV_adj), max(avgtraceMS_adj), max(linsum_adj)])
ax.XTick = [0, xsize/7, (xsize/7)*2, (xsize/7)*3, (xsize/7)*4, (xsize/7)*5, (xsize/7)*6, (xsize/7)*7];
ax.XTickLabel = {'0','1', '2', '3', '4', '5', '6', '7'};
axis([2 xsize-1 y_min y_max])
% print to file (color!!)
clear('fig_filename')
fig_filename=sprintf('F:/Calcium_Imaging_Analysis/tadpoles_byexp/figures/smoothedROIs/paper_examples/adjusted/exp%droi%dsmoothed_averageall.eps', tadpole{1,t}.expnum, r)
% save as eps for illustrator (in color)
saveas(gcf,fig_filename,'psc2');
% save as PNG for easy perusal
clear('fig_filename')
fig_filename=sprintf('F:/Calcium_Imaging_Analysis/tadpoles_byexp/figures/smoothedROIs/paper_examples/adjusted/exp%droi%dsmoothed_averageall.png', tadpole{1,t}.expnum, r)
saveas(gcf,fig_filename,'png');

% exact plots desired for paper:
tracesMS = [2+3, 5+3, 9+3] %
tracesV = [1+3, 5+3, 6+3] %
tracesM = [2+3, 5+3, 9+3] %

% axis dims for individual plots
y_minMS = min([min(filteredMS_adj(tracesMS(1,1),:)), min(filteredMS_adj(tracesMS(1,2),:)), min(filteredMS_adj(tracesMS(1,3),:))])
y_maxMS = max([max(filteredMS_adj(tracesMS(1,1),:)), max(filteredMS_adj(tracesMS(1,2),:)), max(filteredMS_adj(tracesMS(1,3),:))])
y_minM = min([min(filteredM_adj(tracesM(1,1),:)), min(filteredM_adj(tracesM(1,2),:)), min(filteredM_adj(tracesM(1,3),:))])
y_maxM = max([max(filteredM_adj(tracesM(1,1),:)), max(filteredM_adj(tracesM(1,2),:)), max(filteredM_adj(tracesM(1,3),:))])
y_minV = min([min(filteredV_adj(tracesV(1,1),:)), min(filteredV_adj(tracesV(1,2),:)), min(filteredV_adj(tracesV(1,3),:))])
y_maxV = max([max(filteredV_adj(tracesV(1,1),:)), max(filteredV_adj(tracesV(1,2),:)), max(filteredV_adj(tracesV(1,3),:))])
y_min = min([y_minMS, y_minM, y_minV])
y_max = max([y_maxMS, y_maxM, y_maxV])
xsize = size(filteredMS_adj,2);

% multisensory (purple)
figure;
hold on
for k = 1:length(tracesMS)
    plot(filteredMS_adj(tracesMS(k),:), 'color', [0.5 0 0.5])
end
title(sprintf('Exp %d ROI %d smoothed multisensory', tadpole{1,t}.expnum, r))
xlabel('time (s)')
ylabel('\DeltaF/F_{0}')
ax=gca;

%y_min = min([min(filteredMS_adj(tracesMS(1,1),:)), min(filteredMS_adj(tracesMS(1,2),:)), min(filteredMS_adj(tracesMS(1,3),:))])
%y_max = max([max(filteredMS_adj(tracesMS(1,1),:)), max(filteredMS_adj(tracesMS(1,2),:)), max(filteredMS_adj(tracesMS(1,3),:))])
ax.XTick = [0, xsize/7, (xsize/7)*2, (xsize/7)*3, (xsize/7)*4, (xsize/7)*5, (xsize/7)*6, (xsize/7)*7];
ax.XTickLabel = {'0','1', '2', '3', '4', '5', '6', '7'};
hold off
axis([2 xsize-1 y_min y_max])

fig_filename=sprintf('F:/Calcium_Imaging_Analysis/tadpoles_byexp/figures/smoothedROIs/paper_examples/adjusted/exp%droi%dsmoothed_multiselected.eps', tadpole{1,t}.expnum, r)
saveas(gcf,fig_filename,'psc2');

% visual (blue)
figure;
hold on
for k = 1:length(tracesV)
    plot(filteredV_adj(tracesV(k),:), 'color', 'b')
end
title(sprintf('Exp %d ROI %d smoothed visual', tadpole{1,t}.expnum, r))
xlabel('time (s)')
ylabel('\DeltaF/F_{0}')
ax=gca;
xsize = size(filteredMS_adj,2);
%y_min = min([min(filteredMS_adj(tracesMS(1,1),:)), min(filteredMS_adj(tracesMS(1,2),:)), min(filteredMS_adj(tracesMS(1,3),:))])
%y_max = max([max(filteredMS_adj(tracesMS(1,1),:)), max(filteredMS_adj(tracesMS(1,2),:)), max(filteredMS_adj(tracesMS(1,3),:))])
ax.XTick = [0, xsize/7, (xsize/7)*2, (xsize/7)*3, (xsize/7)*4, (xsize/7)*5, (xsize/7)*6, (xsize/7)*7];
ax.XTickLabel = {'0','1', '2', '3', '4', '5', '6', '7'};
hold off
axis([2 xsize-1 y_min y_max])
fig_filename=sprintf('F:/Calcium_Imaging_Analysis/tadpoles_byexp/figures/smoothedROIs/paper_examples/adjusted/exp%droi%dsmoothed_visualselected.eps', tadpole{1,t}.expnum, r)
saveas(gcf,fig_filename,'psc2');

% mechanosensory (red)
figure;
hold on
for k = 1:length(tracesM)
    plot(filteredM_adj(tracesV(k),:), 'color', 'r')
end
title(sprintf('Exp %d ROI %d smoothed mech', tadpole{1,t}.expnum, r))
xlabel('time (s)')
ylabel('\DeltaF/F_{0}')
ax=gca;
ax.XTick = [0, xsize/7, (xsize/7)*2, (xsize/7)*3, (xsize/7)*4, (xsize/7)*5, (xsize/7)*6, (xsize/7)*7];
ax.XTickLabel = {'0','1', '2', '3', '4', '5', '6', '7'};
hold off
axis([2 xsize-1 y_min y_max])
fig_filename=sprintf('F:/Calcium_Imaging_Analysis/tadpoles_byexp/figures/smoothedROIs/paper_examples/adjusted/exp%droi%dsmoothed_mechselected.eps', tadpole{1,t}.expnum, r)
saveas(gcf,fig_filename,'psc2');