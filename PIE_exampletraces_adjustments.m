% % test code for PIE paper examples
% 
% smallest_val = min([min(avgtraceM), min(avgtraceV)])
% avgtraceM_adj = avgtraceM - smallest_val
% avgtraceV_adj = avgtraceV - smallest_val
% avgtraceMS_adj = avgtraceMS - smallest_val
% linsum_adj = avgtraceM_adj + avgtraceV_adj
% 
% % plot all on a graph together
% figure;
% hold on
% plot(avgtraceMS_adj, 'color', [0.5 0 0.5])
% plot(avgtraceV_adj, 'color', 'b')
% plot(avgtraceM_adj, 'color', 'r')
% plot(linsum_adj, 'color', 'k')
% hold off
% title(sprintf('mean responses exp %d roi %d', tadpole{1,t}.expnum, r))
% xlabel('time (s)')
% ylabel('\DeltaF/F_{0}')
% ax=gca;
% xsize = length(tadpole{1,t}.filtered{1,1});
% ax.XTick = [0, xsize/7, (xsize/7)*2, (xsize/7)*3, (xsize/7)*4, (xsize/7)*5, (xsize/7)*6, (xsize/7)*7];
% ax.XTickLabel = {'0','1', '2', '3', '4', '5', '6', '7'};

% figure out slope of decay from dye bleaching using no stim trials
%locate no stim trials
for i = 1:size(stimlist,2)
    filteredNS(i,:) = tadpole{1,t}.filtered{r, stimlist(4,i)};
end
avgtraceNS = mean(filteredNS);
figure;
plot(filteredNS')

% fit a single polynomial (e.g. a line) to each filtered dataset
X_vals = 1:1:160;
for i = 1:size(filteredNS, 1) % over each presentation
    [p(i,:), S(i), mu(i,:)] = polyfit(X_vals, filteredNS(i, :), 1)
    %p(i) = polyfit(X_vals, filteredNS(i, :), 1)
end

avg_slope = mean(p([2:5,7,8],1))
stddev_slope = std(p(:,1))

for i = 1:size(filteredNS, 1)
    figure; plot(filteredNS(i,:))
    title(sprintf('trace %d', i))
end

add_on = -avg_slope * X_vals *0.03

% correct slope values by add_on
avgtraceM_adj = avgtraceM + add_on
avgtraceV_adj = avgtraceV + add_on
avgtraceMS_adj = avgtraceMS + add_on
linsum_adj = avgtraceM_adj + avgtraceV_adj
% plot
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
ax.XTick = [0, xsize/7, (xsize/7)*2, (xsize/7)*3, (xsize/7)*4, (xsize/7)*5, (xsize/7)*6, (xsize/7)*7];
ax.XTickLabel = {'0','1', '2', '3', '4', '5', '6', '7'};


