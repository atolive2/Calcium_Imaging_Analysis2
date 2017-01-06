%% PIE paper example traces. 

%% Small uni/ large multi examples

%% exp 5 roi 51 - responds to all, with multi being large comparatively. 
% exp 5 is in tadpole{1,7}

% plot all traces
figure;
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
% linear sum of unisensory
linsum = avgtraceM + avgtraceV;

% plot all on a graph together
figure;
hold on
plot(avgtraceMS, 'color', [0.5 0 0.5])
plot(avgtraceV, 'color', 'b')
plot(avgtraceM, 'color', 'r')
plot(linsum, 'color', 'k')
hold off
title(sprintf('mean responses exp %d roi %d', tadpole{1,t}.expnum, r))
xlabel('time (s)')
ylabel('\DeltaF/F_{0}')
ax=gca;
xsize = length(tadpole{1,t}.filtered{1,1});
ax.XTick = [0, xsize/7, (xsize/7)*2, (xsize/7)*3, (xsize/7)*4, (xsize/7)*5, (xsize/7)*6, (xsize/7)*7];
ax.XTickLabel = {'0','1', '2', '3', '4', '5', '6', '7'};
axis([0 xsize 0.95 1.05])
% print to file (color!!)
clear('fig_filename')
fig_filename=sprintf('F:/Calcium_Imaging_Analysis/tadpoles_byexp/figures/smoothedROIs/paper_examples/exp%droi%dsmoothed_averageall.eps', tadpole{1,t}.expnum, r)
%print -depsc fig_filename

% save as EPS for illustrator

saveas(gcf,fig_filename,'psc2');
% save as PNG to preserve colors
fig_filename=sprintf('F:/Calcium_Imaging_Analysis/tadpoles_byexp/figures/smoothedROIs/paper_examples/exp%droi%dsmoothed_averageall.png', tadpole{1,t}.expnum, r)
saveas(gcf,fig_filename,'png');

% plot all traces
for i = 1:size(tadpole{1,t}.filtered, 2) % over all trials
    hold on
    if tadpole{1,t}.stimorder(1, i) == 1
        plot(tadpole{1,t}.filtered{r,i}, 'color', [0.5 0 0.5])
    elseif tadpole{1,t}.stimorder(1, i) == 2
        plot(tadpole{1,t}.filtered{r,i}, 'b')
    elseif tadpole{1,t}.stimorder(1, i) == 3
        plot(tadpole{1,t}.filtered{r,i}, 'r')
    elseif tadpole{1,t}.stimorder(1, i) == 4
        plot(tadpole{1,t}.filtered{r,i}, 'color', [0.5 0.5 0.5])
    end
end
title(sprintf('Exp %d ROI %d smoothed all', 5, r))
xlabel('time (s)')
ylabel('\DeltaF/F_{0}')
ax=gca;
xsize = length(tadpole{1,t}.filtered{1,1});
ax.XTick = [0, xsize/7, (xsize/7)*2, (xsize/7)*3, (xsize/7)*4, (xsize/7)*5, (xsize/7)*6, (xsize/7)*7];
ax.XTickLabel = {'0','1', '2', '3', '4', '5', '6', '7'};
hold off

% plot only multisensory traces
figure;
hold on
for i = 1:size(stimlist, 2) % over all trials
    plot(tadpole{1,t}.filtered{r,stimlist(1,i)}, 'color', color_list(i,:))
end
title(sprintf('Exp %d ROI %d smoothed multisensory', tadpole{1,t}.expnum, r))
xlabel('time (s)')
ylabel('\DeltaF/F_{0}')
ax=gca;
xsize = length(tadpole{1,t}.filtered{1,1});
ax.XTick = [0, xsize/7, (xsize/7)*2, (xsize/7)*3, (xsize/7)*4, (xsize/7)*5, (xsize/7)*6, (xsize/7)*7];
ax.XTickLabel = {'0','1', '2', '3', '4', '5', '6', '7'};
hold off
fig_filename=sprintf('F:/Calcium_Imaging_Analysis/tadpoles_byexp/figures/smoothedROIs/paper_examples/exp%droi%dsmoothed_multiall.png', tadpole{1,t}.expnum, r)
saveas(gcf,fig_filename,'png');
%close;
%clear('fig_filename')

% plot only visual traces
figure;
hold on
for i = 1:size(stimlist, 2) % over all trials
    plot(tadpole{1,t}.filtered{r,stimlist(2,i)}, 'color', color_list(i,:))
end
title(sprintf('Exp %d ROI %d smoothed visual', tadpole{1,t}.expnum, r))
xlabel('time (s)')
ylabel('\DeltaF/F_{0}')
ax=gca;
xsize = length(tadpole{1,t}.filtered{1,1});
ax.XTick = [0, xsize/7, (xsize/7)*2, (xsize/7)*3, (xsize/7)*4, (xsize/7)*5, (xsize/7)*6, (xsize/7)*7];
ax.XTickLabel = {'0','1', '2', '3', '4', '5', '6', '7'};
hold off
fig_filename=sprintf('F:/Calcium_Imaging_Analysis/tadpoles_byexp/figures/smoothedROIs/paper_examples/exp%droi%dsmoothed_visualall.png', tadpole{1,t}.expnum, r)
saveas(gcf,fig_filename,'png');
%close;
%clear('fig_filename')

% plot only mechanosensory traces
figure;
hold on
for i = 1:size(stimlist, 2) % over all trials
    plot(tadpole{1,t}.filtered{r,stimlist(3,i)}, 'color', color_list(i,:))
end
title(sprintf('Exp %d ROI %d smoothed mechanosensory', tadpole{1,t}.expnum, r))
xlabel('time (s)')
ylabel('\DeltaF/F_{0}')
ax=gca;
xsize = length(tadpole{1,t}.filtered{1,1});
ax.XTick = [0, xsize/7, (xsize/7)*2, (xsize/7)*3, (xsize/7)*4, (xsize/7)*5, (xsize/7)*6, (xsize/7)*7];
ax.XTickLabel = {'0','1', '2', '3', '4', '5', '6', '7'};
hold off
fig_filename=sprintf('F:/Calcium_Imaging_Analysis/tadpoles_byexp/figures/smoothedROIs/paper_examples/exp%droi%dsmoothed_mechall.png', tadpole{1,t}.expnum, r)
saveas(gcf,fig_filename,'png');
%close;
%clear('fig_filename')

axis([0 xsize -0.08 0.2])

% exact plots desired for paper:
tracesMS = [1, 12, 4] %(Y, C, R)
tracesV = [4, 6, 7] %(R, blue, black)
tracesM = [1, 5, 7] %(Y, gr, black)

% multisensory (purple)
figure;
hold on
for k = 1:length(tracesMS)
    plot(tadpole{1,t}.filtered{r,stimlist(1,tracesMS(k))}, 'color', [0.5 0 0.5])
end
title(sprintf('Exp %d ROI %d smoothed multisensory', tadpole{1,t}.expnum, r))
xlabel('time (s)')
ylabel('\DeltaF/F_{0}')
ax=gca;
xsize = length(tadpole{1,t}.filtered{1,1});
ax.XTick = [0, xsize/7, (xsize/7)*2, (xsize/7)*3, (xsize/7)*4, (xsize/7)*5, (xsize/7)*6, (xsize/7)*7];
ax.XTickLabel = {'0','1', '2', '3', '4', '5', '6', '7'};
hold off
axis([2 xsize -0.02 0.2])

fig_filename=sprintf('F:/Calcium_Imaging_Analysis/tadpoles_byexp/figures/smoothedROIs/paper_examples/exp%droi%dsmoothed_multiselected.eps', tadpole{1,t}.expnum, r)
saveas(gcf,fig_filename,'eps');

% visual (blue)
figure;
hold on
for k = 1:length(tracesV)
    plot(tadpole{1,t}.filtered{r,stimlist(2,tracesV(k))}, 'color', 'b')
end
title(sprintf('Exp %d ROI %d smoothed visual', tadpole{1,t}.expnum, r))
xlabel('time (s)')
ylabel('\DeltaF/F_{0}')
ax=gca;
xsize = length(tadpole{1,t}.filtered{1,1});
ax.XTick = [0, xsize/7, (xsize/7)*2, (xsize/7)*3, (xsize/7)*4, (xsize/7)*5, (xsize/7)*6, (xsize/7)*7];
ax.XTickLabel = {'0','1', '2', '3', '4', '5', '6', '7'};
hold off
axis([2 xsize -0.08 0.2])
fig_filename=sprintf('F:/Calcium_Imaging_Analysis/tadpoles_byexp/figures/smoothedROIs/paper_examples/exp%droi%dsmoothed_visualselected.eps', tadpole{1,t}.expnum, r)
saveas(gcf,fig_filename,'eps');

% mechanosensory (red)
figure;
hold on
for k = 1:length(tracesM)
    plot(tadpole{1,t}.filtered{r,stimlist(3,tracesM(k))}, 'color', 'r')
end
title(sprintf('Exp %d ROI %d smoothed mech', tadpole{1,t}.expnum, r))
xlabel('time (s)')
ylabel('\DeltaF/F_{0}')
ax=gca;
xsize = length(tadpole{1,t}.filtered{1,1});
ax.XTick = [0, xsize/7, (xsize/7)*2, (xsize/7)*3, (xsize/7)*4, (xsize/7)*5, (xsize/7)*6, (xsize/7)*7];
ax.XTickLabel = {'0','1', '2', '3', '4', '5', '6', '7'};
hold off
axis([2 xsize -0.08 0.2])
fig_filename=sprintf('F:/Calcium_Imaging_Analysis/tadpoles_byexp/figures/smoothedROIs/paper_examples/exp%droi%dsmoothed_mechselected.eps', tadpole{1,t}.expnum, r)
saveas(gcf,fig_filename,'eps');


%% Exp 5 roi 151
% exp 5 roi 151 - responds only to multi 
% exp 5 is in tadpole{1,7}

% plot all traces

t = 7;
% identify trials that equal each stim type
clear('stimorder_uniques', 'stimlist')
stimorder_uniques = unique(tadpole{1,t}.stimorder)
for s = 1:length(stimorder_uniques)
    stimlist(s,:) = find(tadpole{1,t}.stimorder == stimorder_uniques(s))
end
color_list = [1 1 0; 1 0 1; 0 1 1; 1 0 0; 0 1 0; 0 0 1; 0 0 0; 0.5 0 0.5; 0.5 0.5 0.5; 0.5 0.5 0; 0 0.5 0.5; 0 0 0.5]     
r = 151;

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
% linear sum of unisensory
linsum = avgtraceM + avgtraceV;

% plot all on a graph together
figure;
hold on
plot(avgtraceMS, 'color', [0.5 0 0.5])
plot(avgtraceV, 'color', 'b')
plot(avgtraceM, 'color', 'r')
plot(linsum, 'color', 'k')
hold off
title(sprintf('mean responses exp %d roi %d', tadpole{1,t}.expnum, r))
xlabel('time (s)')
ylabel('\DeltaF/F_{0}')
ax=gca;
xsize = length(tadpole{1,t}.filtered{1,1});
ax.XTick = [0, xsize/7, (xsize/7)*2, (xsize/7)*3, (xsize/7)*4, (xsize/7)*5, (xsize/7)*6, (xsize/7)*7];
ax.XTickLabel = {'0','1', '2', '3', '4', '5', '6', '7'};
axis([2 xsize-2 -0.04 0.08])
% save as EPS for illustrator
fig_filename=sprintf('F:/Calcium_Imaging_Analysis/tadpoles_byexp/figures/smoothedROIs/paper_examples/exp%droi%dsmoothed_averageall.eps', tadpole{1,t}.expnum, r)
saveas(gcf,fig_filename,'psc2');
% save as PNG to preserve colors
fig_filename=sprintf('F:/Calcium_Imaging_Analysis/tadpoles_byexp/figures/smoothedROIs/paper_examples/exp%droi%dsmoothed_averageall.png', tadpole{1,t}.expnum, r)
saveas(gcf,fig_filename,'png');

% plot all
figure;
for i = 1:size(tadpole{1,t}.filtered, 2) % over all trials
    hold on
    if tadpole{1,t}.stimorder(1, i) == 1
        plot(tadpole{1,t}.filtered{r,i}, 'color', [0.5 0 0.5])
    elseif tadpole{1,t}.stimorder(1, i) == 2
        plot(tadpole{1,t}.filtered{r,i}, 'b')
    elseif tadpole{1,t}.stimorder(1, i) == 3
        plot(tadpole{1,t}.filtered{r,i}, 'r')
    elseif tadpole{1,t}.stimorder(1, i) == 4
        plot(tadpole{1,t}.filtered{r,i}, 'color', [0.5 0.5 0.5])
    end
end
title(sprintf('Exp %d ROI %d smoothed all', 5, r))
xlabel('time (s)')
ylabel('\DeltaF/F_{0}')
ax=gca;
xsize = length(tadpole{1,t}.filtered{1,1});
ax.XTick = [0, xsize/7, (xsize/7)*2, (xsize/7)*3, (xsize/7)*4, (xsize/7)*5, (xsize/7)*6, (xsize/7)*7];
ax.XTickLabel = {'0','1', '2', '3', '4', '5', '6', '7'};
hold off

% plot only multisensory traces
figure;
hold on
for i = 1:size(stimlist, 2) % over all trials
    plot(tadpole{1,t}.filtered{r,stimlist(1,i)}, 'color', color_list(i,:))
end
title(sprintf('Exp %d ROI %d smoothed multisensory', tadpole{1,t}.expnum, r))
xlabel('time (s)')
ylabel('\DeltaF/F_{0}')
ax=gca;
xsize = length(tadpole{1,t}.filtered{1,1});
ax.XTick = [0, xsize/7, (xsize/7)*2, (xsize/7)*3, (xsize/7)*4, (xsize/7)*5, (xsize/7)*6, (xsize/7)*7];
ax.XTickLabel = {'0','1', '2', '3', '4', '5', '6', '7'};
hold off
fig_filename=sprintf('F:/Calcium_Imaging_Analysis/tadpoles_byexp/figures/smoothedROIs/paper_examples/exp%droi%dsmoothed_multiall.png', tadpole{1,t}.expnum, r)
saveas(gcf,fig_filename,'png');
%close;
%clear('fig_filename')

% plot only visual traces
figure;
hold on
for i = 1:size(stimlist, 2) % over all trials
    plot(tadpole{1,t}.filtered{r,stimlist(2,i)}, 'color', color_list(i,:))
end
title(sprintf('Exp %d ROI %d smoothed visual', tadpole{1,t}.expnum, r))
xlabel('time (s)')
ylabel('\DeltaF/F_{0}')
ax=gca;
xsize = length(tadpole{1,t}.filtered{1,1});
ax.XTick = [0, xsize/7, (xsize/7)*2, (xsize/7)*3, (xsize/7)*4, (xsize/7)*5, (xsize/7)*6, (xsize/7)*7];
ax.XTickLabel = {'0','1', '2', '3', '4', '5', '6', '7'};
hold off
fig_filename=sprintf('F:/Calcium_Imaging_Analysis/tadpoles_byexp/figures/smoothedROIs/paper_examples/exp%droi%dsmoothed_visualall.png', tadpole{1,t}.expnum, r)
saveas(gcf,fig_filename,'png');
%close;
%clear('fig_filename')

% plot only mechanosensory traces
figure;
hold on
for i = 1:size(stimlist, 2) % over all trials
    plot(tadpole{1,t}.filtered{r,stimlist(3,i)}, 'color', color_list(i,:))
end
title(sprintf('Exp %d ROI %d smoothed mechanosensory', tadpole{1,t}.expnum, r))
xlabel('time (s)')
ylabel('\DeltaF/F_{0}')
ax=gca;
xsize = length(tadpole{1,t}.filtered{1,1});
ax.XTick = [0, xsize/7, (xsize/7)*2, (xsize/7)*3, (xsize/7)*4, (xsize/7)*5, (xsize/7)*6, (xsize/7)*7];
ax.XTickLabel = {'0','1', '2', '3', '4', '5', '6', '7'};
hold off
fig_filename=sprintf('F:/Calcium_Imaging_Analysis/tadpoles_byexp/figures/smoothedROIs/paper_examples/exp%droi%dsmoothed_mechall.png', tadpole{1,t}.expnum, r)
saveas(gcf,fig_filename,'png');
%close;
%clear('fig_filename')

axis([0 xsize -0.08 0.2])

% exact plots desired for paper:
tracesMS = [1, 4, 6] %(Y, C, R)
tracesV = [2, 3, 5] %(R, blue, black)
tracesM = [1, 2, 4] %(Y, gr, black)

% multisensory (purple)
figure;
hold on
for k = 1:length(tracesMS)
    plot(tadpole{1,t}.filtered{r,stimlist(1,tracesMS(k))}, 'color', [0.5 0 0.5])
end
title(sprintf('Exp %d ROI %d smoothed multisensory', tadpole{1,t}.expnum, r))
xlabel('time (s)')
ylabel('\DeltaF/F_{0}')
ax=gca;
xsize = length(tadpole{1,t}.filtered{1,1});
ax.XTick = [0, xsize/7, (xsize/7)*2, (xsize/7)*3, (xsize/7)*4, (xsize/7)*5, (xsize/7)*6, (xsize/7)*7];
ax.XTickLabel = {'0','1', '2', '3', '4', '5', '6', '7'};
hold off
axis([2 xsize -0.05 0.35])

fig_filename=sprintf('F:/Calcium_Imaging_Analysis/tadpoles_byexp/figures/smoothedROIs/paper_examples/exp%droi%dsmoothed_multiselected.eps', tadpole{1,t}.expnum, r)
saveas(gcf,fig_filename,'eps');

% visual (blue)
figure;
hold on
for k = 1:length(tracesV)
    plot(tadpole{1,t}.filtered{r,stimlist(2,tracesV(k))}, 'color', 'b')
end
title(sprintf('Exp %d ROI %d smoothed visual', tadpole{1,t}.expnum, r))
xlabel('time (s)')
ylabel('\DeltaF/F_{0}')
ax=gca;
xsize = length(tadpole{1,t}.filtered{1,1});
ax.XTick = [0, xsize/7, (xsize/7)*2, (xsize/7)*3, (xsize/7)*4, (xsize/7)*5, (xsize/7)*6, (xsize/7)*7];
ax.XTickLabel = {'0','1', '2', '3', '4', '5', '6', '7'};
hold off
axis([2 xsize -0.08 0.2])
fig_filename=sprintf('F:/Calcium_Imaging_Analysis/tadpoles_byexp/figures/smoothedROIs/paper_examples/exp%droi%dsmoothed_visualselected.eps', tadpole{1,t}.expnum, r)
saveas(gcf,fig_filename,'eps');

% mechanosensory (red)
figure;
hold on
for k = 1:length(tracesM)
    plot(tadpole{1,t}.filtered{r,stimlist(3,tracesM(k))}, 'color', 'r')
end
title(sprintf('Exp %d ROI %d smoothed mech', tadpole{1,t}.expnum, r))
xlabel('time (s)')
ylabel('\DeltaF/F_{0}')
ax=gca;
xsize = length(tadpole{1,t}.filtered{1,1});
ax.XTick = [0, xsize/7, (xsize/7)*2, (xsize/7)*3, (xsize/7)*4, (xsize/7)*5, (xsize/7)*6, (xsize/7)*7];
ax.XTickLabel = {'0','1', '2', '3', '4', '5', '6', '7'};
hold off
axis([2 xsize -0.05 0.35])
fig_filename=sprintf('F:/Calcium_Imaging_Analysis/tadpoles_byexp/figures/smoothedROIs/paper_examples/exp%droi%dsmoothed_mechselected.eps', tadpole{1,t}.expnum, r)
saveas(gcf,fig_filename,'eps');

%% Exp 6 roi 57
% exp 6 roi 57 - responds only to multi and V
% exp 6 is in tadpole{1,8}

% plot all traces
t = 8;
% identify trials that equal each stim type
clear('stimorder_uniques', 'stimlist')
stimorder_uniques = unique(tadpole{1,t}.stimorder)
for s = 1:length(stimorder_uniques)
    stimlist(s,:) = find(tadpole{1,t}.stimorder == stimorder_uniques(s))
end
color_list = [1 1 0; 1 0 1; 0 1 1; 1 0 0; 0 1 0; 0 0 1; 0 0 0; 0.5 0 0.5; 0.5 0.5 0.5; 0.5 0.5 0; 0 0.5 0.5; 0 0 0.5]     
r = 57;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% averages for each modality
% generate matrix of traces I need
clear('filteredMS', 'filteredM', 'filteredV', 'avgtraceMS', 'avgtraceM', 'avgtraceV')
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
% linear sum of unisensory
linsum = avgtraceM + avgtraceV;

% plot all on a graph together
figure;
hold on
plot(avgtraceMS, 'color', [0.5 0 0.5])
plot(avgtraceV, 'color', 'b')
plot(avgtraceM, 'color', 'r')
plot(linsum, 'color', 'k')
hold off
title(sprintf('mean responses exp %d roi %d', tadpole{1,t}.expnum, r))
xlabel('time (s)')
ylabel('\DeltaF/F_{0}')
ax=gca;
xsize = length(tadpole{1,t}.filtered{1,1});
ax.XTick = [0, xsize/7, (xsize/7)*2, (xsize/7)*3, (xsize/7)*4, (xsize/7)*5, (xsize/7)*6, (xsize/7)*7];
ax.XTickLabel = {'0','1', '2', '3', '4', '5', '6', '7'};
axis([0 xsize-1 -0.04 0.1])
% save as EPS for illustrator
fig_filename=sprintf('F:/Calcium_Imaging_Analysis/tadpoles_byexp/figures/smoothedROIs/paper_examples/exp%droi%dsmoothed_averageall.eps', tadpole{1,t}.expnum, r)
saveas(gcf,fig_filename,'psc2');
% save as PNG to preserve colors
fig_filename=sprintf('F:/Calcium_Imaging_Analysis/tadpoles_byexp/figures/smoothedROIs/paper_examples/exp%droi%dsmoothed_averageall.png', tadpole{1,t}.expnum, r)
saveas(gcf,fig_filename,'png');

% plot all traces
figure;
for i = 1:size(tadpole{1,t}.filtered, 2) % over all trials
    hold on
    if tadpole{1,t}.stimorder(1, i) == 1
        plot(tadpole{1,t}.filtered{r,i}, 'color', [0.5 0 0.5])
    elseif tadpole{1,t}.stimorder(1, i) == 2
        plot(tadpole{1,t}.filtered{r,i}, 'b')
    elseif tadpole{1,t}.stimorder(1, i) == 3
        plot(tadpole{1,t}.filtered{r,i}, 'r')
    elseif tadpole{1,t}.stimorder(1, i) == 4
        plot(tadpole{1,t}.filtered{r,i}, 'color', [0.5 0.5 0.5])
    end
end
title(sprintf('Exp %d ROI %d smoothed all', 5, r))
xlabel('time (s)')
ylabel('\DeltaF/F_{0}')
ax=gca;
xsize = length(tadpole{1,t}.filtered{1,80});
ax.XTick = [0, xsize/7, (xsize/7)*2, (xsize/7)*3, (xsize/7)*4, (xsize/7)*5, (xsize/7)*6, (xsize/7)*7];
ax.XTickLabel = {'0','1', '2', '3', '4', '5', '6', '7'};
hold off
axis([0 xsize -0.1 0.5])

% plot only multisensory traces
figure;
hold on
for i = 13:size(stimlist, 2) % over all trials
    plot(tadpole{1,t}.filtered{r,stimlist(1,i)}, 'color', color_list((i-12),:))
end
title(sprintf('Exp %d ROI %d smoothed multisensory', tadpole{1,t}.expnum, r))
xlabel('time (s)')
ylabel('\DeltaF/F_{0}')
ax=gca;
xsize = length(tadpole{1,t}.filtered{1,80});
ax.XTick = [0, xsize/7, (xsize/7)*2, (xsize/7)*3, (xsize/7)*4, (xsize/7)*5, (xsize/7)*6, (xsize/7)*7];
ax.XTickLabel = {'0','1', '2', '3', '4', '5', '6', '7'};
hold off
fig_filename=sprintf('F:/Calcium_Imaging_Analysis/tadpoles_byexp/figures/smoothedROIs/paper_examples/exp%droi%dsmoothed_multiall.png', tadpole{1,t}.expnum, r)
saveas(gcf,fig_filename,'png');
%close;
%clear('fig_filename')

% plot only visual traces
figure;
hold on
for i = 13:size(stimlist, 2) % over all trials
    plot(tadpole{1,t}.filtered{r,stimlist(2,i)}, 'color', color_list((i-12),:))
end
title(sprintf('Exp %d ROI %d smoothed visual', tadpole{1,t}.expnum, r))
xlabel('time (s)')
ylabel('\DeltaF/F_{0}')
ax=gca;
xsize = length(tadpole{1,t}.filtered{1,80});
ax.XTick = [0, xsize/7, (xsize/7)*2, (xsize/7)*3, (xsize/7)*4, (xsize/7)*5, (xsize/7)*6, (xsize/7)*7];
ax.XTickLabel = {'0','1', '2', '3', '4', '5', '6', '7'};
hold off
fig_filename=sprintf('F:/Calcium_Imaging_Analysis/tadpoles_byexp/figures/smoothedROIs/paper_examples/exp%droi%dsmoothed_visualall.png', tadpole{1,t}.expnum, r)
saveas(gcf,fig_filename,'png');
%close;
%clear('fig_filename')

% plot only mechanosensory traces
figure;
hold on
for i = 13:size(stimlist, 2) % over all trials
    plot(tadpole{1,t}.filtered{r,stimlist(3,i)}, 'color', color_list(i-12,:))
end
title(sprintf('Exp %d ROI %d smoothed mechanosensory', tadpole{1,t}.expnum, r))
xlabel('time (s)')
ylabel('\DeltaF/F_{0}')
ax=gca;
xsize = length(tadpole{1,t}.filtered{1,80});
ax.XTick = [0, xsize/7, (xsize/7)*2, (xsize/7)*3, (xsize/7)*4, (xsize/7)*5, (xsize/7)*6, (xsize/7)*7];
ax.XTickLabel = {'0','1', '2', '3', '4', '5', '6', '7'};
hold off
fig_filename=sprintf('F:/Calcium_Imaging_Analysis/tadpoles_byexp/figures/smoothedROIs/paper_examples/exp%droi%dsmoothed_mechall.png', tadpole{1,t}.expnum, r)
saveas(gcf,fig_filename,'png');
%close;
%clear('fig_filename')

axis([0 xsize -0.2 0.5])

% exact plots desired for paper:
tracesMS = [6+12, 8+12, 9+12] %
tracesV = [2+12, 8+12, 9+12] %
tracesM = [2+12, 5+12, 8+12] %

% multisensory (purple)
figure;
hold on
for k = 1:length(tracesMS)
    plot(tadpole{1,t}.filtered{r,stimlist(1,tracesMS(k))}, 'color', [0.5 0 0.5])
end
title(sprintf('Exp %d ROI %d smoothed multisensory', tadpole{1,t}.expnum, r))
xlabel('time (s)')
ylabel('\DeltaF/F_{0}')
ax=gca;
xsize = length(tadpole{1,t}.filtered{1,80});
ax.XTick = [0, xsize/7, (xsize/7)*2, (xsize/7)*3, (xsize/7)*4, (xsize/7)*5, (xsize/7)*6, (xsize/7)*7];
ax.XTickLabel = {'0','1', '2', '3', '4', '5', '6', '7'};
hold off
axis([2 xsize-2 -0.2 0.5])

fig_filename=sprintf('F:/Calcium_Imaging_Analysis/tadpoles_byexp/figures/smoothedROIs/paper_examples/exp%droi%dsmoothed_multiselected.eps', tadpole{1,t}.expnum, r)
saveas(gcf,fig_filename,'eps');

% visual (blue)
figure;
hold on
for k = 1:length(tracesV)
    plot(tadpole{1,t}.filtered{r,stimlist(2,tracesV(k))}, 'color', 'b')
end
title(sprintf('Exp %d ROI %d smoothed visual', tadpole{1,t}.expnum, r))
xlabel('time (s)')
ylabel('\DeltaF/F_{0}')
ax=gca;
xsize = length(tadpole{1,t}.filtered{1,80});
ax.XTick = [0, xsize/7, (xsize/7)*2, (xsize/7)*3, (xsize/7)*4, (xsize/7)*5, (xsize/7)*6, (xsize/7)*7];
ax.XTickLabel = {'0','1', '2', '3', '4', '5', '6', '7'};
hold off
axis([2 xsize-2 -0.2 0.5])
fig_filename=sprintf('F:/Calcium_Imaging_Analysis/tadpoles_byexp/figures/smoothedROIs/paper_examples/exp%droi%dsmoothed_visualselected.eps', tadpole{1,t}.expnum, r)
saveas(gcf,fig_filename,'eps');

% mechanosensory (red)
figure;
hold on
for k = 1:length(tracesM)
    plot(tadpole{1,t}.filtered{r,stimlist(3,tracesM(k))}, 'color', 'r')
end
title(sprintf('Exp %d ROI %d smoothed mech', tadpole{1,t}.expnum, r))
xlabel('time (s)')
ylabel('\DeltaF/F_{0}')
ax=gca;
xsize = length(tadpole{1,t}.filtered{1,80});
ax.XTick = [0, xsize/7, (xsize/7)*2, (xsize/7)*3, (xsize/7)*4, (xsize/7)*5, (xsize/7)*6, (xsize/7)*7];
ax.XTickLabel = {'0','1', '2', '3', '4', '5', '6', '7'};
hold off
axis([2 xsize-2 -0.2 0.5])
fig_filename=sprintf('F:/Calcium_Imaging_Analysis/tadpoles_byexp/figures/smoothedROIs/paper_examples/exp%droi%dsmoothed_mechselected.eps', tadpole{1,t}.expnum, r)
saveas(gcf,fig_filename,'eps');

%% large uni/same multi

%% Exp 8 roi 19 (V largest)
% stored in tadpole{1,10}

% plot all traces
figure;
t = 10;
% identify trials that equal each stim type
stimorder_uniques = unique(tadpole{1,t}.stimorder)
clear('stimlist')
for s = 1:length(stimorder_uniques)
    stimlist(s,:) = find(tadpole{1,t}.stimorder == stimorder_uniques(s))
end
color_list = [1 1 0; 1 0 1; 0 1 1; 1 0 0; 0 1 0; 0 0 1; 0 0 0; 0.5 0 0.5; 0.5 0.5 0.5; 0.5 0.5 0; 0 0.5 0.5; 0 0 0.5]     
r =19;

for i = 1:size(tadpole{1,t}.filtered, 2) % over all trials
    hold on
    if tadpole{1,t}.stimorder(1, i) == 1
        plot(tadpole{1,t}.filtered{r,i}, 'color', [0.5 0 0.5])
    elseif tadpole{1,t}.stimorder(1, i) == 2
        plot(tadpole{1,t}.filtered{r,i}, 'b')
    elseif tadpole{1,t}.stimorder(1, i) == 3
        plot(tadpole{1,t}.filtered{r,i}, 'r')
    elseif tadpole{1,t}.stimorder(1, i) == 4
        plot(tadpole{1,t}.filtered{r,i}, 'color', [0.5 0.5 0.5])
    end
end
title(sprintf('Exp %d ROI %d smoothed all', 5, r))
xlabel('time (s)')
ylabel('\DeltaF/F_{0}')
ax=gca;
xsize = length(tadpole{1,t}.filtered{1,40});
ax.XTick = [0, xsize/7, (xsize/7)*2, (xsize/7)*3, (xsize/7)*4, (xsize/7)*5, (xsize/7)*6, (xsize/7)*7];
ax.XTickLabel = {'0','1', '2', '3', '4', '5', '6', '7'};
hold off
axis([0 xsize -0.1 0.5])

% plot only multisensory traces
figure;
hold on
for i = 4:size(stimlist, 2) % over all trials
    plot(tadpole{1,t}.filtered{r,stimlist(1,i)}, 'color', color_list((i-3),:))
end
title(sprintf('Exp %d ROI %d smoothed multisensory', tadpole{1,t}.expnum, r))
xlabel('time (s)')
ylabel('\DeltaF/F_{0}')
ax=gca;
xsize = length(tadpole{1,t}.filtered{1,40});
ax.XTick = [0, xsize/7, (xsize/7)*2, (xsize/7)*3, (xsize/7)*4, (xsize/7)*5, (xsize/7)*6, (xsize/7)*7];
ax.XTickLabel = {'0','1', '2', '3', '4', '5', '6', '7'};
hold off
fig_filename=sprintf('F:/Calcium_Imaging_Analysis/tadpoles_byexp/figures/smoothedROIs/paper_examples/exp%droi%dsmoothed_multiall.png', tadpole{1,t}.expnum, r)
saveas(gcf,fig_filename,'png');
%close;
%clear('fig_filename')

% plot only visual traces
figure;
hold on
for i = 4:size(stimlist, 2) % over all trials
    plot(tadpole{1,t}.filtered{r,stimlist(2,i)}, 'color', color_list((i-3),:))
end
title(sprintf('Exp %d ROI %d smoothed visual', tadpole{1,t}.expnum, r))
xlabel('time (s)')
ylabel('\DeltaF/F_{0}')
ax=gca;
xsize = length(tadpole{1,t}.filtered{1,80});
ax.XTick = [0, xsize/7, (xsize/7)*2, (xsize/7)*3, (xsize/7)*4, (xsize/7)*5, (xsize/7)*6, (xsize/7)*7];
ax.XTickLabel = {'0','1', '2', '3', '4', '5', '6', '7'};
hold off
fig_filename=sprintf('F:/Calcium_Imaging_Analysis/tadpoles_byexp/figures/smoothedROIs/paper_examples/exp%droi%dsmoothed_visualall.png', tadpole{1,t}.expnum, r)
saveas(gcf,fig_filename,'png');
%close;
%clear('fig_filename')

% plot only mechanosensory traces
figure;
hold on
for i = 4:size(stimlist, 2) % over all trials
    plot(tadpole{1,t}.filtered{r,stimlist(3,i)}, 'color', color_list(i-3,:))
end
title(sprintf('Exp %d ROI %d smoothed mechanosensory', tadpole{1,t}.expnum, r))
xlabel('time (s)')
ylabel('\DeltaF/F_{0}')
ax=gca;
xsize = length(tadpole{1,t}.filtered{1,80});
ax.XTick = [0, xsize/7, (xsize/7)*2, (xsize/7)*3, (xsize/7)*4, (xsize/7)*5, (xsize/7)*6, (xsize/7)*7];
ax.XTickLabel = {'0','1', '2', '3', '4', '5', '6', '7'};
hold off
fig_filename=sprintf('F:/Calcium_Imaging_Analysis/tadpoles_byexp/figures/smoothedROIs/paper_examples/exp%droi%dsmoothed_mechall.png', tadpole{1,t}.expnum, r)
saveas(gcf,fig_filename,'png');
%close;
%clear('fig_filename')

axis([0 xsize -0.1 0.3])

% exact plots desired for paper:
tracesMS = [2+3, 5+3, 9+3] %
tracesV = [1+3, 5+3, 6+3] %
tracesM = [2+3, 5+3, 9+3] %

% multisensory (purple)
figure;
hold on
for k = 1:length(tracesMS)
    plot(tadpole{1,t}.filtered{r,stimlist(1,tracesMS(k))}, 'color', [0.5 0 0.5])
end
title(sprintf('Exp %d ROI %d smoothed multisensory', tadpole{1,t}.expnum, r))
xlabel('time (s)')
ylabel('\DeltaF/F_{0}')
ax=gca;
xsize = length(tadpole{1,t}.filtered{1,40});
ax.XTick = [0, xsize/7, (xsize/7)*2, (xsize/7)*3, (xsize/7)*4, (xsize/7)*5, (xsize/7)*6, (xsize/7)*7];
ax.XTickLabel = {'0','1', '2', '3', '4', '5', '6', '7'};
hold off
axis([2 xsize-2 -0.2 0.5])

fig_filename=sprintf('F:/Calcium_Imaging_Analysis/tadpoles_byexp/figures/smoothedROIs/paper_examples/exp%droi%dsmoothed_multiselected.eps', tadpole{1,t}.expnum, r)
saveas(gcf,fig_filename,'eps');

% visual (blue)
figure;
hold on
for k = 1:length(tracesV)
    plot(tadpole{1,t}.filtered{r,stimlist(2,tracesV(k))}, 'color', 'b')
end
title(sprintf('Exp %d ROI %d smoothed visual', tadpole{1,t}.expnum, r))
xlabel('time (s)')
ylabel('\DeltaF/F_{0}')
ax=gca;
xsize = length(tadpole{1,t}.filtered{1,40});
ax.XTick = [0, xsize/7, (xsize/7)*2, (xsize/7)*3, (xsize/7)*4, (xsize/7)*5, (xsize/7)*6, (xsize/7)*7];
ax.XTickLabel = {'0','1', '2', '3', '4', '5', '6', '7'};
hold off
axis([2 xsize-2 -0.2 0.5])
fig_filename=sprintf('F:/Calcium_Imaging_Analysis/tadpoles_byexp/figures/smoothedROIs/paper_examples/exp%droi%dsmoothed_visualselected.eps', tadpole{1,t}.expnum, r)
saveas(gcf,fig_filename,'eps');

% mechanosensory (red)
figure;
hold on
for k = 1:length(tracesM)
    plot(tadpole{1,t}.filtered{r,stimlist(3,tracesM(k))}, 'color', 'r')
end
title(sprintf('Exp %d ROI %d smoothed mech', tadpole{1,t}.expnum, r))
xlabel('time (s)')
ylabel('\DeltaF/F_{0}')
ax=gca;
xsize = length(tadpole{1,t}.filtered{1,40});
ax.XTick = [0, xsize/7, (xsize/7)*2, (xsize/7)*3, (xsize/7)*4, (xsize/7)*5, (xsize/7)*6, (xsize/7)*7];
ax.XTickLabel = {'0','1', '2', '3', '4', '5', '6', '7'};
hold off
axis([2 xsize-2 -0.06 0.3])
fig_filename=sprintf('F:/Calcium_Imaging_Analysis/tadpoles_byexp/figures/smoothedROIs/paper_examples/exp%droi%dsmoothed_mechselected.eps', tadpole{1,t}.expnum, r)
saveas(gcf,fig_filename,'eps');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% averages for each modality
% generate matrix of traces I need
clear('filteredMS', 'filteredM', 'filteredV', 'avgtraceMS', 'avgtraceM', 'avgtraceV')
% multisensory
for i = 4:size(stimlist,2)
    filteredMS(i-3,:) = tadpole{1,t}.filtered{r, stimlist(1,i)};
end
avgtraceMS = mean(filteredMS)
% visual
for i = 4:size(stimlist,2)
    filteredV(i-3,:) = tadpole{1,t}.filtered{r, stimlist(2,i)};
end
avgtraceV = mean(filteredV);
%mechanosensory
for i = 4:size(stimlist,2)
    filteredM(i-3,:) = tadpole{1,t}.filtered{r, stimlist(3,i)};
end
avgtraceM = mean(filteredM);
% linear sum of unisensory
linsum = avgtraceM + avgtraceV;

% plot all on a graph together
figure;
hold on
plot(avgtraceMS, 'color', [0.5 0 0.5])
plot(avgtraceV, 'color', 'b')
plot(avgtraceM, 'color', 'r')
plot(linsum, 'color', 'k')
hold off
title(sprintf('mean responses exp %d roi %d', tadpole{1,t}.expnum, r))
xlabel('time (s)')
ylabel('\DeltaF/F_{0}')
ax=gca;
xsize = length(tadpole{1,t}.filtered{1,40});
ax.XTick = [0, xsize/7, (xsize/7)*2, (xsize/7)*3, (xsize/7)*4, (xsize/7)*5, (xsize/7)*6, (xsize/7)*7];
ax.XTickLabel = {'0','1', '2', '3', '4', '5', '6', '7'};
% save as EPS for illustrator
fig_filename=sprintf('F:/Calcium_Imaging_Analysis/tadpoles_byexp/figures/smoothedROIs/paper_examples/exp%droi%dsmoothed_averageall.eps', tadpole{1,t}.expnum, r)
saveas(gcf,fig_filename,'psc2');
% save as PNG to preserve colors
fig_filename=sprintf('F:/Calcium_Imaging_Analysis/tadpoles_byexp/figures/smoothedROIs/paper_examples/exp%droi%dsmoothed_averageall.png', tadpole{1,t}.expnum, r)
saveas(gcf,fig_filename,'png');

%% Exp 8 roi 22 (very large multisensory)
% stored in tadpole{1,10}

% plot all traces

t = 10;
% identify trials that equal each stim type
stimorder_uniques = unique(tadpole{1,t}.stimorder)
clear('stimlist')
for s = 1:length(stimorder_uniques)
    stimlist(s,:) = find(tadpole{1,t}.stimorder == stimorder_uniques(s))
end
color_list = [1 1 0; 1 0 1; 0 1 1; 1 0 0; 0 1 0; 0 0 1; 0 0 0; 0.5 0 0.5; 0.5 0.5 0.5; 0.5 0.5 0; 0 0.5 0.5; 0 0 0.5]     
r =22;

figure;
for i = 1:size(tadpole{1,t}.filtered, 2) % over all trials
    hold on
    if tadpole{1,t}.stimorder(1, i) == 1
        plot(tadpole{1,t}.filtered{r,i}, 'color', [0.5 0 0.5])
    elseif tadpole{1,t}.stimorder(1, i) == 2
        plot(tadpole{1,t}.filtered{r,i}, 'b')
    elseif tadpole{1,t}.stimorder(1, i) == 3
        plot(tadpole{1,t}.filtered{r,i}, 'r')
    elseif tadpole{1,t}.stimorder(1, i) == 4
        plot(tadpole{1,t}.filtered{r,i}, 'color', [0.5 0.5 0.5])
    end
end
title(sprintf('Exp %d ROI %d smoothed all', 5, r))
xlabel('time (s)')
ylabel('\DeltaF/F_{0}')
ax=gca;
xsize = length(tadpole{1,t}.filtered{1,40});
ax.XTick = [0, xsize/7, (xsize/7)*2, (xsize/7)*3, (xsize/7)*4, (xsize/7)*5, (xsize/7)*6, (xsize/7)*7];
ax.XTickLabel = {'0','1', '2', '3', '4', '5', '6', '7'};
hold off
axis([0 xsize -0.1 0.5])

% plot only multisensory traces
figure;
hold on
for i = 4:size(stimlist, 2) % over all trials
    plot(tadpole{1,t}.filtered{r,stimlist(1,i)}, 'color', color_list((i-3),:))
end
title(sprintf('Exp %d ROI %d smoothed multisensory', tadpole{1,t}.expnum, r))
xlabel('time (s)')
ylabel('\DeltaF/F_{0}')
ax=gca;
xsize = length(tadpole{1,t}.filtered{1,40});
ax.XTick = [0, xsize/7, (xsize/7)*2, (xsize/7)*3, (xsize/7)*4, (xsize/7)*5, (xsize/7)*6, (xsize/7)*7];
ax.XTickLabel = {'0','1', '2', '3', '4', '5', '6', '7'};
hold off
fig_filename=sprintf('F:/Calcium_Imaging_Analysis/tadpoles_byexp/figures/smoothedROIs/paper_examples/exp%droi%dsmoothed_multiall.png', tadpole{1,t}.expnum, r)
saveas(gcf,fig_filename,'png');
%close;
%clear('fig_filename')

% plot only visual traces
figure;
hold on
for i = 4:size(stimlist, 2) % over all trials
    plot(tadpole{1,t}.filtered{r,stimlist(2,i)}, 'color', color_list((i-3),:))
end
title(sprintf('Exp %d ROI %d smoothed visual', tadpole{1,t}.expnum, r))
xlabel('time (s)')
ylabel('\DeltaF/F_{0}')
ax=gca;
xsize = length(tadpole{1,t}.filtered{1,40});
ax.XTick = [0, xsize/7, (xsize/7)*2, (xsize/7)*3, (xsize/7)*4, (xsize/7)*5, (xsize/7)*6, (xsize/7)*7];
ax.XTickLabel = {'0','1', '2', '3', '4', '5', '6', '7'};
hold off
fig_filename=sprintf('F:/Calcium_Imaging_Analysis/tadpoles_byexp/figures/smoothedROIs/paper_examples/exp%droi%dsmoothed_visualall.png', tadpole{1,t}.expnum, r)
saveas(gcf,fig_filename,'png');
%close;
%clear('fig_filename')

% plot only mechanosensory traces
figure;
hold on
for i = 4:size(stimlist, 2) % over all trials
    plot(tadpole{1,t}.filtered{r,stimlist(3,i)}, 'color', color_list(i-3,:))
end
title(sprintf('Exp %d ROI %d smoothed mechanosensory', tadpole{1,t}.expnum, r))
xlabel('time (s)')
ylabel('\DeltaF/F_{0}')
ax=gca;
xsize = length(tadpole{1,t}.filtered{1,40});
ax.XTick = [0, xsize/7, (xsize/7)*2, (xsize/7)*3, (xsize/7)*4, (xsize/7)*5, (xsize/7)*6, (xsize/7)*7];
ax.XTickLabel = {'0','1', '2', '3', '4', '5', '6', '7'};
hold off
fig_filename=sprintf('F:/Calcium_Imaging_Analysis/tadpoles_byexp/figures/smoothedROIs/paper_examples/exp%droi%dsmoothed_mechall.png', tadpole{1,t}.expnum, r)
saveas(gcf,fig_filename,'png');
%close;
%clear('fig_filename')

axis([0 xsize -0.1 0.5])

% exact plots desired for paper:
tracesMS = [5+3, 3+3, 8+3] %
tracesV = [4+3, 6+3, 5+3] %
tracesM = [5+3, 6+3, 8+3] %

% multisensory (purple)
figure;
hold on
for k = 1:length(tracesMS)
    plot(tadpole{1,t}.filtered{r,stimlist(1,tracesMS(k))}, 'color', [0.5 0 0.5])
end
title(sprintf('Exp %d ROI %d smoothed multisensory', tadpole{1,t}.expnum, r))
xlabel('time (s)')
ylabel('\DeltaF/F_{0}')
ax=gca;
xsize = length(tadpole{1,t}.filtered{1,40});
ax.XTick = [0, xsize/7, (xsize/7)*2, (xsize/7)*3, (xsize/7)*4, (xsize/7)*5, (xsize/7)*6, (xsize/7)*7];
ax.XTickLabel = {'0','1', '2', '3', '4', '5', '6', '7'};
hold off
axis([2 xsize-2 -0.1 0.5])

fig_filename=sprintf('F:/Calcium_Imaging_Analysis/tadpoles_byexp/figures/smoothedROIs/paper_examples/exp%droi%dsmoothed_multiselected.eps', tadpole{1,t}.expnum, r)
saveas(gcf,fig_filename,'eps');

% visual (blue)
figure;
hold on
for k = 1:length(tracesV)
    plot(tadpole{1,t}.filtered{r,stimlist(2,tracesV(k))}, 'color', 'b')
end
title(sprintf('Exp %d ROI %d smoothed visual', tadpole{1,t}.expnum, r))
xlabel('time (s)')
ylabel('\DeltaF/F_{0}')
ax=gca;
xsize = length(tadpole{1,t}.filtered{1,40});
ax.XTick = [0, xsize/7, (xsize/7)*2, (xsize/7)*3, (xsize/7)*4, (xsize/7)*5, (xsize/7)*6, (xsize/7)*7];
ax.XTickLabel = {'0','1', '2', '3', '4', '5', '6', '7'};
hold off
axis([2 xsize-2 -0.1 0.5])
fig_filename=sprintf('F:/Calcium_Imaging_Analysis/tadpoles_byexp/figures/smoothedROIs/paper_examples/exp%droi%dsmoothed_visualselected.eps', tadpole{1,t}.expnum, r)
saveas(gcf,fig_filename,'eps');

% mechanosensory (red)
figure;
hold on
for k = 1:length(tracesM)
    plot(tadpole{1,t}.filtered{r,stimlist(3,tracesM(k))}, 'color', 'r')
end
title(sprintf('Exp %d ROI %d smoothed mech', tadpole{1,t}.expnum, r))
xlabel('time (s)')
ylabel('\DeltaF/F_{0}')
ax=gca;
xsize = length(tadpole{1,t}.filtered{1,40});
ax.XTick = [0, xsize/7, (xsize/7)*2, (xsize/7)*3, (xsize/7)*4, (xsize/7)*5, (xsize/7)*6, (xsize/7)*7];
ax.XTickLabel = {'0','1', '2', '3', '4', '5', '6', '7'};
hold off
axis([2 xsize-2 -0.1 0.5])
fig_filename=sprintf('F:/Calcium_Imaging_Analysis/tadpoles_byexp/figures/smoothedROIs/paper_examples/exp%droi%dsmoothed_mechselected.eps', tadpole{1,t}.expnum, r)
saveas(gcf,fig_filename,'eps');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% averages for each modality
% generate matrix of traces I need
clear('filteredMS', 'filteredM', 'filteredV', 'avgtraceMS', 'avgtraceM', 'avgtraceV')
% multisensory
for i = 4:size(stimlist,2)
    filteredMS(i-3,:) = tadpole{1,t}.filtered{r, stimlist(1,i)};
end
avgtraceMS = mean(filteredMS)
% visual
for i = 4:size(stimlist,2)
    filteredV(i-3,:) = tadpole{1,t}.filtered{r, stimlist(2,i)};
end
avgtraceV = mean(filteredV);
%mechanosensory
for i = 4:size(stimlist,2)
    filteredM(i-3,:) = tadpole{1,t}.filtered{r, stimlist(3,i)};
end
avgtraceM = mean(filteredM);
% linear sum of unisensory
linsum = avgtraceM + avgtraceV;

% plot all on a graph together
figure;
hold on
plot(avgtraceMS, 'color', [0.5 0 0.5])
plot(avgtraceV, 'color', 'b')
plot(avgtraceM, 'color', 'r')
plot(linsum, 'color', 'k')
hold off
title(sprintf('mean responses exp %d roi %d', tadpole{1,t}.expnum, r))
xlabel('time (s)')
ylabel('\DeltaF/F_{0}')
ax=gca;
xsize = length(tadpole{1,t}.filtered{1,40});
ax.XTick = [0, xsize/7, (xsize/7)*2, (xsize/7)*3, (xsize/7)*4, (xsize/7)*5, (xsize/7)*6, (xsize/7)*7];
ax.XTickLabel = {'0','1', '2', '3', '4', '5', '6', '7'};
% save as EPS for illustrator
fig_filename=sprintf('F:/Calcium_Imaging_Analysis/tadpoles_byexp/figures/smoothedROIs/paper_examples/exp%droi%dsmoothed_averageall.eps', tadpole{1,t}.expnum, r)
saveas(gcf,fig_filename,'psc2');
% save as PNG to preserve colors
fig_filename=sprintf('F:/Calcium_Imaging_Analysis/tadpoles_byexp/figures/smoothedROIs/paper_examples/exp%droi%dsmoothed_averageall.png', tadpole{1,t}.expnum, r)
saveas(gcf,fig_filename,'png');

%% Exp 6, roi 45 (V/M/MS all the same size (roughly)

% exp 6 is in tadpole{1,8}
t = 8;
r = 45;


% plot all traces
figure;
for i = 1:size(tadpole{1,t}.filtered, 2) % over all trials
    hold on
    if tadpole{1,t}.stimorder(1, i) == 1
        plot(tadpole{1,t}.filtered{r,i}, 'color', [0.5 0 0.5])
    elseif tadpole{1,t}.stimorder(1, i) == 2
        plot(tadpole{1,t}.filtered{r,i}, 'b')
    elseif tadpole{1,t}.stimorder(1, i) == 3
        plot(tadpole{1,t}.filtered{r,i}, 'r')
    elseif tadpole{1,t}.stimorder(1, i) == 4
        plot(tadpole{1,t}.filtered{r,i}, 'color', [0.5 0.5 0.5])
    end
end
title(sprintf('Exp %d ROI %d smoothed all', 5, r))
xlabel('time (s)')
ylabel('\DeltaF/F_{0}')
ax=gca;
xsize = length(tadpole{1,t}.filtered{1,80});
ax.XTick = [0, xsize/7, (xsize/7)*2, (xsize/7)*3, (xsize/7)*4, (xsize/7)*5, (xsize/7)*6, (xsize/7)*7];
ax.XTickLabel = {'0','1', '2', '3', '4', '5', '6', '7'};
hold off
axis([0 xsize -0.1 0.5])

% identify trials that equal each stim type
clear('stimorder_uniques')
stimorder_uniques = unique(tadpole{1,t}.stimorder)
clear('stimlist')
for s = 1:length(stimorder_uniques)
    stimlist(s,:) = find(tadpole{1,t}.stimorder == stimorder_uniques(s))
end
color_list = [1 1 0; 1 0 1; 0 1 1; 1 0 0; 0 1 0; 0 0 1; 0 0 0; 0.5 0 0.5; 0.5 0.5 0.5; 0.5 0.5 0; 0 0.5 0.5; 0 0 0.5]     

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% averages for each modality
% generate matrix of traces I need
clear('filteredMS', 'filteredM', 'filteredV', 'avgtraceMS', 'avgtraceM', 'avgtraceV')
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
% linear sum of unisensory
linsum = avgtraceM + avgtraceV;

% plot all on a graph together
figure;
hold on
plot(avgtraceMS, 'color', [0.5 0 0.5])
plot(avgtraceV, 'color', 'b')
plot(avgtraceM, 'color', 'r')
plot(linsum, 'color', 'k')
hold off
title(sprintf('mean responses exp %d roi %d', tadpole{1,t}.expnum, r))
xlabel('time (s)')
ylabel('\DeltaF/F_{0}')
ax=gca;
xsize = length(tadpole{1,t}.filtered{1,1});
ax.XTick = [0, xsize/7, (xsize/7)*2, (xsize/7)*3, (xsize/7)*4, (xsize/7)*5, (xsize/7)*6, (xsize/7)*7];
ax.XTickLabel = {'0','1', '2', '3', '4', '5', '6', '7'};
axis([2 xsize-1 -0.05 0.25])
% save as EPS for illustrator
fig_filename=sprintf('F:/Calcium_Imaging_Analysis/tadpoles_byexp/figures/smoothedROIs/paper_examples/exp%droi%dsmoothed_averageall.epsc', tadpole{1,t}.expnum, r)
saveas(gcf,fig_filename,'epsc');

print -depsc fig_filename

% save as PNG to preserve colors
fig_filename=sprintf('F:/Calcium_Imaging_Analysis/tadpoles_byexp/figures/smoothedROIs/paper_examples/exp%droi%dsmoothed_averageall.png', tadpole{1,t}.expnum, r)
saveas(gcf,fig_filename,'png');


% plot only multisensory traces
figure;
hold on
for i = 1:6 % over all trials
    plot(tadpole{1,t}.filtered{r,stimlist(1,i)}, 'color', color_list((i),:))
end
title(sprintf('Exp %d ROI %d smoothed multisensory', tadpole{1,t}.expnum, r))
xlabel('time (s)')
ylabel('\DeltaF/F_{0}')
ax=gca;
xsize = length(tadpole{1,t}.filtered{1,80});
ax.XTick = [0, xsize/7, (xsize/7)*2, (xsize/7)*3, (xsize/7)*4, (xsize/7)*5, (xsize/7)*6, (xsize/7)*7];
ax.XTickLabel = {'0','1', '2', '3', '4', '5', '6', '7'};
hold off
fig_filename=sprintf('F:/Calcium_Imaging_Analysis/tadpoles_byexp/figures/smoothedROIs/paper_examples/exp%droi%dsmoothed_multiall.png', tadpole{1,t}.expnum, r)
saveas(gcf,fig_filename,'png');
%close;
%clear('fig_filename')

% plot only visual traces
figure;
hold on
for i = 1:6 % over all trials
    plot(tadpole{1,t}.filtered{r,stimlist(2,i)}, 'color', color_list((i),:))
end
title(sprintf('Exp %d ROI %d smoothed visual', tadpole{1,t}.expnum, r))
xlabel('time (s)')
ylabel('\DeltaF/F_{0}')
ax=gca;
xsize = length(tadpole{1,t}.filtered{1,80});
ax.XTick = [0, xsize/7, (xsize/7)*2, (xsize/7)*3, (xsize/7)*4, (xsize/7)*5, (xsize/7)*6, (xsize/7)*7];
ax.XTickLabel = {'0','1', '2', '3', '4', '5', '6', '7'};
hold off
fig_filename=sprintf('F:/Calcium_Imaging_Analysis/tadpoles_byexp/figures/smoothedROIs/paper_examples/exp%droi%dsmoothed_visualall.png', tadpole{1,t}.expnum, r)
saveas(gcf,fig_filename,'png');
%close;
%clear('fig_filename')

% plot only mechanosensory traces
figure;
hold on
for i = 1:6 % over all trials
    plot(tadpole{1,t}.filtered{r,stimlist(3,i)}, 'color', color_list(i,:))
end
title(sprintf('Exp %d ROI %d smoothed mechanosensory', tadpole{1,t}.expnum, r))
xlabel('time (s)')
ylabel('\DeltaF/F_{0}')
ax=gca;
xsize = length(tadpole{1,t}.filtered{1,80});
ax.XTick = [0, xsize/7, (xsize/7)*2, (xsize/7)*3, (xsize/7)*4, (xsize/7)*5, (xsize/7)*6, (xsize/7)*7];
ax.XTickLabel = {'0','1', '2', '3', '4', '5', '6', '7'};
hold off
fig_filename=sprintf('F:/Calcium_Imaging_Analysis/tadpoles_byexp/figures/smoothedROIs/paper_examples/exp%droi%dsmoothed_mechall.png', tadpole{1,t}.expnum, r)
saveas(gcf,fig_filename,'png');
%close;
%clear('fig_filename')

axis([0 xsize -0.2 0.5])

% exact plots desired for paper:
tracesMS = [2, 4, 6] %
tracesV = [2, 4, 5] %
tracesM = [1, 5, 6] %

% multisensory (purple)
figure;
hold on
for k = 1:length(tracesMS)
    plot(tadpole{1,t}.filtered{r,stimlist(1,tracesMS(k))}, 'color', [0.5 0 0.5])
end
title(sprintf('Exp %d ROI %d smoothed multisensory', tadpole{1,t}.expnum, r))
xlabel('time (s)')
ylabel('\DeltaF/F_{0}')
ax=gca;
xsize = length(tadpole{1,t}.filtered{1,5});
ax.XTick = [0, xsize/7, (xsize/7)*2, (xsize/7)*3, (xsize/7)*4, (xsize/7)*5, (xsize/7)*6, (xsize/7)*7];
ax.XTickLabel = {'0','1', '2', '3', '4', '5', '6', '7'};
hold off
axis([2 xsize-2 -0.2 0.5])

fig_filename=sprintf('F:/Calcium_Imaging_Analysis/tadpoles_byexp/figures/smoothedROIs/paper_examples/exp%droi%dsmoothed_multiselected.eps', tadpole{1,t}.expnum, r)
saveas(gcf,fig_filename,'eps');

% visual (blue)
figure;
hold on
for k = 1:length(tracesV)
    plot(tadpole{1,t}.filtered{r,stimlist(2,tracesV(k))}, 'color', 'b')
end
title(sprintf('Exp %d ROI %d smoothed visual', tadpole{1,t}.expnum, r))
xlabel('time (s)')
ylabel('\DeltaF/F_{0}')
ax=gca;
xsize = length(tadpole{1,t}.filtered{1,5});
ax.XTick = [0, xsize/7, (xsize/7)*2, (xsize/7)*3, (xsize/7)*4, (xsize/7)*5, (xsize/7)*6, (xsize/7)*7];
ax.XTickLabel = {'0','1', '2', '3', '4', '5', '6', '7'};
hold off
axis([2 xsize-2 -0.2 0.5])
fig_filename=sprintf('F:/Calcium_Imaging_Analysis/tadpoles_byexp/figures/smoothedROIs/paper_examples/exp%droi%dsmoothed_visualselected.eps', tadpole{1,t}.expnum, r)
saveas(gcf,fig_filename,'eps');

% mechanosensory (red)
figure;
hold on
for k = 1:length(tracesM)
    plot(tadpole{1,t}.filtered{r,stimlist(3,tracesM(k))}, 'color', 'r')
end
title(sprintf('Exp %d ROI %d smoothed mech', tadpole{1,t}.expnum, r))
xlabel('time (s)')
ylabel('\DeltaF/F_{0}')
ax=gca;
xsize = length(tadpole{1,t}.filtered{1,4});
ax.XTick = [0, xsize/7, (xsize/7)*2, (xsize/7)*3, (xsize/7)*4, (xsize/7)*5, (xsize/7)*6, (xsize/7)*7];
ax.XTickLabel = {'0','1', '2', '3', '4', '5', '6', '7'};
hold off
axis([2 xsize-2 -0.1 0.35])
fig_filename=sprintf('F:/Calcium_Imaging_Analysis/tadpoles_byexp/figures/smoothedROIs/paper_examples/exp%droi%dsmoothed_mechselected.eps', tadpole{1,t}.expnum, r)
saveas(gcf,fig_filename,'eps');