%% Dissertation figure 5: xcorr basics

%% 5A - example correlated cell
% use example from AAAS 2018 poster
% tad 34, roi 2 vs 7, 8, 50, 63, 64, 67, multisensory trials 4-15

% what t is tad 34?
for t = 1:length(allData)
    tads(t) = allData{1,t}.expnum
end
t = find(tads == 34)
% t is 20

% plot ROI 2 vs 7, 8, 50, 63

% use allData_old to get dffo_multi field (see tad34_examplePlots for
% original code)
%rois = [2 7 8 50 63];
rois = [2 63 8 7 50]
% try again
figure;
k=0
hold on
for r = 1:length(rois)
    plot(allData_old{1,20}.dff0_multi(rois(r), :)+k, 'Color', H(2, rois(r),:), 'LineWidth', 2);
    k = k-0.5
end
hold off

xlim([636 2385])
%ylim([-1 1])
ax = gca;
%ax.Visible = 'off'
%title(sprintf('Tad %d xcorr maxR with ROI %d', t, r))
ylabel('\DeltaF/F_{0}')
xlabel('time (frames)')
fig = gcf;
fig.PaperUnits = 'inches'
fig.PaperPosition = [0 0 8 3]
saveas(gcf, 'tad34_roi2vs63,8,7,50', 'png')
saveas(gcf, 'tad34_roi2vs63,8,7,50', 'epsc2')

%% Example xcorr plots

% B - show 0 lag R ~= maxR
% tad 43 (t=25), st46, V.

% C - show 46 more stripy
% tads 48 (t=28) (46) and 3 (t=2) (49)

% D - show 49 more square
% tads 44 (t=26) 46, 31 (t=17) 49

%% 5E: quantify stripe/square/gradient



%% Notes

for r = 1:length(tadcluster{1,t}.resp_ROIs)
    if r == roi
        continue
    else
        figure;
    hold on
    
    plot(tadcluster{1,t}.dff0_multi(roi,:), 'Color', 'k', 'LineWidth', 2)
    plot(tadcluster{1,t}.dff0_multi(r,:)+0.3, 'Color', H(roi,r,:), 'LineWidth', 2)

    hold off
    px = round(pdist2(tadcluster{1,t}.ROIcenters(roi,:), tadcluster{1,t}.ROIcenters(r,:), 'euclidean'), 0)
    xlim([636 2385])
    ylim([-1 1])
    ax = gca;
    ax.Visible = 'off'
%     title(sprintf('Tad %d ROI %d vs ROI %d (%d px apart)', t, roi, r, px), 'fontsize', 30)
%     ylabel('\DeltaF/F_{0}', 'fontsize', 30)
%     xlabel('time (frames)', 'fontsize', 30)
%     set(gca, 'fontsize', 30)
    fig = gcf;
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 8 3];
    fig_filename = sprintf('Tad %d ROI %d vs ROI %d (%d px apart) no axes', tadcluster{1,t}.expnum, roi, r, px)
    saveas(gcf, fig_filename, 'png')