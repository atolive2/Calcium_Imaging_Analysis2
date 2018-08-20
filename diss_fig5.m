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

% Calculate average R val for maxR and 0lagR for MS, V, M
for t = 1:length(allData)
    if isfield(allData{1,t}, 'respROIdff0_maxR_sq_MS')
        tmp = find(allData{1,t}.respROIdff0_maxR_sq_MS > 0 & allData{1,t}.respROIdff0_maxR_sq_MS < 1);
        avg_Rmax(t,1) = mean(allData{1,t}.respROIdff0_maxR_sq_MS(tmp));
        clear('tmp')
        tmp = find(allData{1,t}.respROIdff0_maxR_sq_V > 0 & allData{1,t}.respROIdff0_maxR_sq_V < 1);
        avg_Rmax(t,2) = mean(allData{1,t}.respROIdff0_maxR_sq_V(tmp));
        clear('tmp')
        tmp = find(allData{1,t}.respROIdff0_maxR_sq_M > 0 & allData{1,t}.respROIdff0_maxR_sq_M < 1);
        avg_Rmax(t,3) = mean(allData{1,t}.respROIdff0_maxR_sq_M(tmp));
        clear('tmp')
    else 
        avg_Rmax(t,1) = NaN
        avg_Rmax(t,2) = NaN
        avg_Rmax(t,3) = NaN
    end
end

% Get mean value by stage (in plotting order)
%46 first
for m = 1:size(avg_Rmax,2)
    mean_avg_Rmax(m) = nanmean(avg_Rmax(s46_tads,m));
end
%49 second
for m = 1:size(avg_Rmax,2)
    mean_avg_Rmax(m+3) = nanmean(avg_Rmax(s49_tads,m));
end

%sd
%46
for m = 1:size(avg_Rmax,2)
    mean_avg_Rmax(2,m) = nanstd(avg_Rmax(s46_tads,m));
end
%49 second
for m = 1:size(avg_Rmax,2)
    mean_avg_Rmax(2,m+3) = nanstd(avg_Rmax(s49_tads,m));
end

% Make plot (same as proportion HC ROIs)
st46_X = ones(1, length(s46_tads))
st49_X = 2*ones(1, length(s49_tads))
m=15
st = 1

figure;
hold on
plot(st46_X, avg_Rmax(s46_tads, 1), 'oy', 'MarkerSize', m)%46 MS 
plot(st46_X+0.3, avg_Rmax(s46_tads, 2), 'or', 'MarkerSize', m) %46 V
plot(st46_X+0.6, avg_Rmax(s46_tads, 3), 'ob', 'MarkerSize', m) %46 MS
plot(st49_X, avg_Rmax(s49_tads, 1), 'oy', 'MarkerSize', m) %49 MS 
plot(st49_X+0.3, avg_Rmax(s49_tads, 2), 'or', 'MarkerSize', m) %49 V
plot(st49_X+0.6, avg_Rmax(s49_tads, 3), 'ob', 'MarkerSize', m) %49 MS

plot([(st-0.1) (st+0.1)], [mean_avg_Rmax(1,1), mean_avg_Rmax(1,1)], 'k', 'LineWidth', 6) %46 MS
plot([(st-0.1+0.3) (st+0.1+0.3)], [mean_avg_Rmax(1,2), mean_avg_Rmax(1,2)], 'k', 'LineWidth', 6) % 46 V
plot([(st-0.1+0.6) (st+0.1+0.6)], [mean_avg_Rmax(3), mean_avg_Rmax(3)], 'k', 'LineWidth', 6) % 46 V

plot([(st-0.1+1) (st+0.1+1)], [mean_avg_Rmax(1,4), mean_avg_Rmax(1,4)], 'k', 'LineWidth', 6) %49 MS
plot([(st-0.1+1.3) (st+0.1+1.3)], [mean_avg_Rmax(1,5), mean_avg_Rmax(1,5)], 'k', 'LineWidth', 6) % 49 V
plot([(st-0.1+1.6) (st+0.1+1.6)], [mean_avg_Rmax(1,6), mean_avg_Rmax(1,6)], 'k', 'LineWidth', 6) % 49 V

hold off 
xlim([0.6 3])
ax = gca;
ax.XTick = [1 1.3 1.6 2 2.3 2.6]
ax.XTickLabel = {'MS', 'V', 'M', 'MS', 'V', 'M'}
ylabel('mean corr')
set(gca, 'FontSize', 30)
saveas(gcf, 'average xcorr by mod and st', 'png')
saveas(gcf, 'average xcorr by mod and st', 'epsc2')

% collect all maxR avgs into a matrix for Carlos
for t = 1:length(allData)
    key(t) = allData{1,t}.stage
end
s46_tads = find(key == 46)
s49_tads = find(key == 49)

% stage 46
avg_Rmax_46 = avg_Rmax(s46_tads, :)
avg_Rmax_49 = avg_Rmax(s49_tads, :)

%% stats on avg correlation
% do an ANOVA using anovan (unbalanced design)
% data vector (remove tads 5 and 16 b/c not enough data)

data = [avg_Rmax(:,1); avg_Rmax(:,2); avg_Rmax(:,3)];
data(:,2) = [ones(size(avg_Rmax,1),1); 2*ones(size(avg_Rmax,1),1); 3*ones(size(avg_Rmax,1),1)];
data(:,3) = [key'; key';, key'];
% remove 5 and 16
is = ~isnan(data(:,1));
data_cl = data(is,:)
% remove discoboxed tad 19
isnot1 = find(data_cl(:,3) ~= 1);
data_cl = data_cl(isnot1,:)

% Run anovan
[p, tbl, stats] = anovan(data_cl(:,1), {data_cl(:,2), data_cl(:,3)})
% mod p = 0.54, stage p = 0.15

%% Stats on distributions of correlations

% Do a Kruskal-Wallis test on each tad for MS vs V vs M.
% KW appears to be the only nonparametric test to see if 3+ distributions are not pulled
% from the same population. 
% cite: https://www.mathworks.com/help/stats/kruskalwallis.html#btv4oqy-1-displayopt

for t = 1:length(allData)
    %get data
    if isfield(allData{1,t}, 'respROIdff0_maxR_sq_MS')
        tmp = allData{1,t}.respROIdff0_maxR_sq_MS(find(allData{1,t}.respROIdff0_maxR_sq_MS > 0 & allData{1,t}.respROIdff0_maxR_sq_MS < .9999));
        tmp1 = allData{1,t}.respROIdff0_maxR_sq_V(find(allData{1,t}.respROIdff0_maxR_sq_V > 0 & allData{1,t}.respROIdff0_maxR_sq_V < .9999));
        tmp2 = allData{1,t}.respROIdff0_maxR_sq_M(find(allData{1,t}.respROIdff0_maxR_sq_M > 0 & allData{1,t}.respROIdff0_maxR_sq_M < .9999)); 
        grp = [ones(length(tmp),1); 2*ones(length(tmp1),1); 3*ones(length(tmp2),1)];
        [maxRdist_P(t), maxRdist_tbl{t}, maxRdist_stats{t}] = kruskalwallis([tmp; tmp1; tmp2], grp, 'displayopt', 'off')
        allData{1,t}.maxR_fordist = [[tmp; tmp1; tmp2], grp]
        clear('tmp', 'tmp1', 'tmp2')
    end
end

%% Plot the P values by stage
P46 = maxRdist_P(st46)
P49 = maxRdist_P(st49)
Pmed(1) = nanmedian(maxRdist_P(st46))
Pmed(2) = nanmedian(maxRdist_P(st49))
st46_X = ones(length(P46),1)
figure;
hold on 
scatter(st46_X, P46, 80, 'g', 'LineWidth', 1)
scatter(st49_X - 0.5, P49, 80, 'm', 'LineWidth', 1)
plot([0.9 1.1], [Pmed(1), Pmed(1)], 'k', 'LineWidth', 2)
plot([1.4 1.6], [Pmed(2), Pmed(2)], 'k', 'LineWidth', 2)
hold off
xlim([0.5 2])
%ylim([0 1.1])
ax = gca;
ax.XTick = [1 1.5];
ax.XTickLabel = [46 49];
%title('Proportion of Responding ROIs by Exp')
xlabel('Stage')
ylabel('P_{KW}')
set(gca, 'FontSize', 30)
saveas(gcf, 'P_KW of by mod maxRval by stage all', 'epsc2')

% zoom in on sig Ps
figure;
hold on 
scatter(st46_X, P46, 80, 'g', 'LineWidth', 1)
scatter(st49_X - 0.5, P49, 80, 'm', 'LineWidth', 1)
plot([0.9 1.1], [Pmed(1), Pmed(1)], 'k', 'LineWidth', 2)
plot([1.4 1.6], [Pmed(2), Pmed(2)], 'k', 'LineWidth', 2)
hold off
xlim([0.5 2])
ylim([0 0.05])
ax = gca;
ax.XTick = [1 1.5];
ax.XTickLabel = [46 49];
%title('Proportion of Responding ROIs by Exp')
xlabel('Stage')
ylabel('P_{KW}')
set(gca, 'FontSize', 30)
saveas(gcf, 'P_KW of by mod maxRval by stage zoomin', 'epsc2')



%% Make histograms for each tad, overlay MS, V, M

for t = 1:length(allData)
    if isfield(allData{1,t}, 'maxR_fordist')
        if ~isempty(allData{1,t}.maxR_fordist)
    MS = find(allData{1,t}.maxR_fordist(:,2) == 1);
    V = find(allData{1,t}.maxR_fordist(:,2) == 2);
    M = find(allData{1,t}.maxR_fordist(:,2) == 3);
    meds(1) = nanmedian(allData{1,t}.maxR_fordist(MS, 1));
    meds(2) = nanmedian(allData{1,t}.maxR_fordist(V, 1));
    meds(3) = nanmedian(allData{1,t}.maxR_fordist(M, 1));
    figure;
    hold on
    histogram(allData{1,t}.maxR_fordist(MS, 1), 'FaceColor', [255/256,215/256,0], 'EdgeColor', [255/256,215/256,0], 'NumBins', 40, 'BinLimits', [0 1], 'Normalization', 'probability');
    histogram(allData{1,t}.maxR_fordist(V, 1), 'FaceColor', 'r', 'EdgeColor', 'r', 'NumBins', 40, 'BinLimits', [0 1], 'Normalization', 'probability');
    histogram(allData{1,t}.maxR_fordist(M, 1), 'FaceColor', 'b', 'EdgeColor', 'b', 'NumBins', 40, 'BinLimits', [0 1], 'Normalization', 'probability');
    plot([meds(1) meds(1)], ylim, 'Color', [255/256,215/256,0], 'LineWidth', 2)
    plot([meds(2) meds(2)], ylim, 'r', 'LineWidth', 2)
    plot([meds(3) meds(3)], ylim, 'Color', [ 0.5843 0.8157 0.9882], 'LineWidth', 2) %light blue so visible
    hold off
    xlabel('max R val')
    ylabel('ROI count')
    fig_filename = sprintf('maxR all hist stg %d tad %d (t=%d)', allData{1,t}.stage, allData{1,t}.expnum, t);
    set(gca, 'FontSize', 30)
    saveas(gcf, fig_filename, 'png')
    saveas(gcf, fig_filename, 'epsc2')    
    clear('MS', 'V', 'M', 'meds')
        end
    end
end

%% Make a histogram of all data for each mod by stage
% collect all data for each stage and mod into vectors
MS_all46 = [];
V_all46 = [];
M_all46 = [];
MS_all49 = [];
V_all49= [];
M_all49 = [];
for t = 1:length(allData)
    if isfield(allData{1,t}, 'maxR_fordist')
        if ~isempty(allData{1,t}.maxR_fordist)
            MS = find(allData{1,t}.maxR_fordist(:,2) == 1);
            V = find(allData{1,t}.maxR_fordist(:,2) == 2);
            M = find(allData{1,t}.maxR_fordist(:,2) == 3);
            if allData{1,t}.stage == 46
                MS_all46 = [MS_all46; allData{1,t}.maxR_fordist(MS, 1)];
                V_all46 = [V_all46; allData{1,t}.maxR_fordist(V, 1)];
                M_all46 = [M_all46; allData{1,t}.maxR_fordist(M, 1)];
            elseif allData{1,t}.stage == 49
                MS_all49 = [MS_all49; allData{1,t}.maxR_fordist(MS, 1)];
                V_all49 = [V_all49; allData{1,t}.maxR_fordist(V, 1)];
                M_all49 = [M_all49; allData{1,t}.maxR_fordist(M, 1)];
            end
        end
    end
end

% calculate median
meds_all(1) = nanmedian(MS_all46)
meds_all(2) = nanmedian(V_all46)
meds_all(3) = nanmedian(M_all46)
meds_all(4) = nanmedian(MS_all49)
meds_all(5) = nanmedian(V_all49)
meds_all(6) = nanmedian(M_all49)

% stats (KStest between stages for each modality)
[h(1), p(1)] = kstest2(MS_all46, MS_all49) % p = 0
[h(2), p(2)] = kstest2(V_all46, V_all49) % p < 0.001
[h(3), p(3)] = kstest2(M_all46, M_all49) % p < 0.001

% Histogram each modality, both stages
% MS
figure;
hold on
histogram(MS_all46, 'FaceColor', 'g', 'EdgeColor', 'g', 'NumBins', 40, 'BinLimits', [0 1], 'Normalization', 'probability');
histogram(MS_all49, 'FaceColor', 'm', 'EdgeColor', 'm', 'NumBins', 40, 'BinLimits', [0 1], 'Normalization', 'probability');
plot([meds_all(1) meds_all(1)], ylim, 'g', 'LineWidth', 2)
plot([meds_all(4) meds_all(4)], ylim, 'm', 'LineWidth', 2)
hold off
xlabel('max R val')
ylabel('norm ROI count')
fig_filename = 'maxR all hist MS by stage'
set(gca, 'FontSize', 30)
saveas(gcf, fig_filename, 'png')
saveas(gcf, fig_filename, 'epsc2')
    
% V
figure;
hold on
histogram(V_all46, 'FaceColor', 'g', 'EdgeColor', 'g', 'NumBins', 40, 'BinLimits', [0 1], 'Normalization', 'probability');
histogram(V_all49, 'FaceColor', 'm', 'EdgeColor', 'm', 'NumBins', 40, 'BinLimits', [0 1], 'Normalization', 'probability');
plot([meds_all(2) meds_all(2)], ylim, 'g', 'LineWidth', 2)
plot([meds_all(5) meds_all(5)], ylim, 'm', 'LineWidth', 2)
hold off
xlabel('max R val')
ylabel('norm ROI count')
fig_filename = 'maxR all hist V by stage'
set(gca, 'FontSize', 30)
saveas(gcf, fig_filename, 'png')
saveas(gcf, fig_filename, 'epsc2')
    
% M
figure;
hold on
histogram(M_all46, 'FaceColor', 'g', 'EdgeColor', 'g', 'NumBins', 40, 'BinLimits', [0 1], 'Normalization', 'probability');
histogram(M_all49, 'FaceColor', 'm', 'EdgeColor', 'm', 'NumBins', 40, 'BinLimits', [0 1], 'Normalization', 'probability');
plot([meds_all(3) meds_all(3)], ylim, 'g', 'LineWidth', 2)
plot([meds_all(6) meds_all(6)], ylim, 'm', 'LineWidth', 2)
hold off
xlabel('max R val')
ylabel('norm ROI count')
fig_filename = 'maxR all hist M by stage'
set(gca, 'FontSize', 30)
saveas(gcf, fig_filename, 'png')
saveas(gcf, fig_filename, 'epsc2')
    
% Histogram all mods, each stage
%46
figure;
hold on
histogram(MS_all46, 'FaceColor', [255/256,215/256,0], 'EdgeColor', [255/256,215/256,0], 'NumBins', 40, 'BinLimits', [0 1], 'Normalization', 'probability');
histogram(V_all46, 'FaceColor', 'r', 'EdgeColor', 'r', 'NumBins', 40, 'BinLimits', [0 1], 'Normalization', 'probability');
histogram(M_all46, 'FaceColor', 'b', 'EdgeColor', 'b', 'NumBins', 40, 'BinLimits', [0 1], 'Normalization', 'probability');
plot([meds_all(1) meds_all(1)], ylim, 'Color', [255/256,215/256,0], 'LineWidth', 2)
plot([meds_all(2) meds_all(2)], ylim, 'r', 'LineWidth', 2)
plot([meds_all(3) meds_all(3)], ylim, 'Color', [ 0.5843 0.8157 0.9882], 'LineWidth', 2) %light blue so visible
hold off
xlabel('max R val')
ylabel('norm ROI count')
fig_filename = 'all tads hist MS V M st46';
set(gca, 'FontSize', 30)
saveas(gcf, fig_filename, 'png')
saveas(gcf, fig_filename, 'epsc2')   

%49
figure;
hold on
histogram(MS_all49, 'FaceColor', [255/256,215/256,0], 'EdgeColor', [255/256,215/256,0], 'NumBins', 40, 'BinLimits', [0 1], 'Normalization', 'probability');
histogram(V_all49, 'FaceColor', 'r', 'EdgeColor', 'r', 'NumBins', 40, 'BinLimits', [0 1], 'Normalization', 'probability');
histogram(M_all49, 'FaceColor', 'b', 'EdgeColor', 'b', 'NumBins', 40, 'BinLimits', [0 1], 'Normalization', 'probability');
plot([meds_all(4) meds_all(4)], ylim, 'Color', [255/256,215/256,0], 'LineWidth', 2)
plot([meds_all(5) meds_all(5)], ylim, 'r', 'LineWidth', 2)
plot([meds_all(6) meds_all(6)], ylim, 'Color', [ 0.5843 0.8157 0.9882], 'LineWidth', 2) %light blue so visible
hold off
xlabel('max R val')
ylabel('norm ROI count')
fig_filename = 'all tads hist MS V M st46';
set(gca, 'FontSize', 30)
saveas(gcf, fig_filename, 'png')
saveas(gcf, fig_filename, 'epsc2')   

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