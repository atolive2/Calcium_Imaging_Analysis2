%% Figure 1B

counts = [0 1 2 3 4]

% get peak data by spike count
%st 46
tmp_peaks = max([CellsAll46(:).trace]);
allpeaks46 = [CellsAll46(:).earlySpikes; CellsAll46(:).spikeCount; tmp_peaks];
for c = 1:length(counts)
    include = find((allpeaks46(1,:) == 0) & (allpeaks46(2,:) == counts(c)))
    peaks_toplot{c,1} = allpeaks46(3, include)
end
clear('tmp_peaks')
%st 49
for t = 1:length(CellsAll49)
    tmp_peaks(t) = max([CellsAll49(t).trace]);
end
%tmp_peaks = [CellsAll49(:).trace(1:159)];
allpeaks49 = [CellsAll49(:).earlySpikes; CellsAll49(:).spikeCount; tmp_peaks];
for c = 1:length(counts)
    include = find((allpeaks49(1,:) == 0) & (allpeaks49(2,:) == counts(c)))
    peaks_toplot{c,2} = allpeaks49(3, include)
end

% get means and error bars
% run an ANOVA using anovan (for unbalanced n)
% get data into right shape
Y = [];
g1 = [];
g2 = [];
%st 46
for p = 1:size(peaks_toplot,1)
    Y = [Y, peaks_toplot{p, 1}];
    g1 = [g1, 46*ones(1, length(peaks_toplot{p, 1}))];
    g2 = [g2, (p-1)*ones(1, length(peaks_toplot{p, 1}))];
end
%st 49
for p = 1:size(peaks_toplot,1)
    Y = [Y, peaks_toplot{p, 2}];
    g1 = [g1, 49*ones(1, length(peaks_toplot{p, 2}))];
    g2 = [g2, (p-1)*ones(1, length(peaks_toplot{p, 2}))];
end

[p, tbl, stats] = anovan(Y, {g1, g2}, 'model','interaction','varnames',{'stage','spike_ct'})
% sig di
% make error barsff across spike count, not across stage. 
c = multcompare(stats, 'Dimension', [1 2])

for c = 1:size(peaks_toplot,1)
    for d = 1:size(peaks_toplot, 2)
        avg_val(c,d) = mean(peaks_toplot{c,d});
        error_val(c,d) = std(peaks_toplot{c,d}) / sqrt(length(peaks_toplot{c,d}));
    end
end
% resize
avg_vals = [];
err_vals = [];
for c = 1:size(peaks_toplot,1)
    avg_vals = [avg_vals, avg_val(c, :)];
    err_vals = [err_vals, error_val(c,:)];
end


%% Make scatter plot
m=10
x_vals = rand(1,85) *0.1;
err_xvals = [-0.1, 0.1 0.9 1.1 1.9 2.1 2.9 3.1 3.9 4.1];
figure;
hold on
plot((x_vals - 0.1), peaks_toplot{1,1}, 'go', 'MarkerSize', m) %46, 0
plot((x_vals(1:71) + 0.1), peaks_toplot{1,2}, 'mo', 'MarkerSize', m) %49, 0
plot((x_vals(1:60) - 0.1 + 1), peaks_toplot{2,1}, 'go', 'MarkerSize', m) %46, 1
plot((x_vals(1:77) + 0.1 + 1), peaks_toplot{2,2}, 'ms', 'MarkerSize', m) %49, 1
plot((x_vals(1:59) - 0.1 + 2), peaks_toplot{3,1}, 'go', 'MarkerSize', m) %46, 2
plot((x_vals(1:51) + 0.1 + 2), peaks_toplot{3,2}, 'ms', 'MarkerSize', m) %49, 2
plot((x_vals(1:39) - 0.1 + 3), peaks_toplot{4,1}, 'go', 'MarkerSize', m) %46, 3
plot((x_vals(1:24) + 0.1 + 3), peaks_toplot{4,2}, 'ms', 'MarkerSize', m) %49, 3
plot((x_vals(1:13) - 0.1 + 4), peaks_toplot{5,1}, 'go', 'MarkerSize', m) %46, 4
plot((x_vals(1:33) + 0.1 + 4), peaks_toplot{5,2}, 'ms', 'MarkerSize', m) %49, 4
e = errorbar(err_xvals, avg_vals, err_vals, '+', 'MarkerSize', 10, 'MarkerEdgeColor', 'red', 'MarkerFaceColor', 'red') %, 'Linewidth', 2, 'LineColor', 'red')
e.Color = 'red'
hold off

xlabel('Spike Count')
ylabel('Peak \DeltaF/F_{0}')
ax = gca;
ax.XTick = [0 1 2 3 4]
xlim([-0.2 4.2])
set(gca, 'FontSize', 30)
% current folder is F:/Dissertation_images
saveas(gcf, 'peak_vs_spikenum', 'png')
saveas(gcf, 'peak_vs_spikenum', 'epsc2')


%% Updated versions (figures v2)

% re run stats as 2 1 way ANOVAs 
% get data into right shape
%stage 46
Y46 = [];
g146 = [];
g246 = [];
for p = 1:size(peaks_toplot,1)
    Y46 = [Y46, peaks_toplot{p, 1}];
    g146 = [g146, 46*ones(1, length(peaks_toplot{p, 1}))];
    g246 = [g246, (p-1)*ones(1, length(peaks_toplot{p, 1}))];
end
[p46, tbl46, stats46] = anova1(Y46, g246) %, 'model','interaction')
[c46, m46] = multcompare(stats46)

%st 49
Y49 = [];
g149 = [];
g249 = [];
for p = 1:size(peaks_toplot,1)
    Y49 = [Y49, peaks_toplot{p, 2}];
    g149 = [g149, 49*ones(1, length(peaks_toplot{p, 2}))];
    g249 = [g249, (p-1)*ones(1, length(peaks_toplot{p, 2}))];
end

[p49, tbl49, stats49] = anova1(Y49, g249) %, 'model','interaction')
[c49, m49] = multcompare(stats49)

%all data together
Y_all = [];
g1_all = [];
g2_all = [];
for p = 1:size(peaks_toplot,1)
    Y_all = [Y_all, peaks_toplot{p, 2}];
    g1_all = [g1_all, 49*ones(1, length(peaks_toplot{p, 2}))];
    g2_all = [g2_all, (p-1)*ones(1, length(peaks_toplot{p, 2}))];
end

[p_all, tb_all, stats_all] = anova1(Y_all, g2_all) %, 'model','interaction')
[c_all, m_all] = multcompare(stats_all)

% ANOVAs are all signnificant - so there is a relationship over spike
% count. 

%% Create 0 spikes mean and conf int (to show 0.1 threshold)
counter = 1
for i = 1:length(CellsAll49)
    if CellsAll49(i).spikeCount == 0
        all_0spikes(:,counter) = CellsAll49(i).trace(1:159);
        counter = counter + 1
    end
end

avg_0spikes = mean(all_0spikes,2);
std_0spikes = std(all_0spikes');
confint_0spikes = avg_0spikes + -

CI = mean(x)+- t * (s / square(n)) 


figure;
hold on
plot(all_0spikes, 'Color', [0.5 0.5 0.5])
plot(avg_0spikes, 'k', 'LineWidth', 3)

%% 0 spikes mean and conf int using 46s
counter = 1
for i = 1:length(CellsAll46)
    if CellsAll46(i).spikeCount == 0
        all_0spikes46(:,counter) = CellsAll46(i).trace(1:159);
        counter = counter + 1
    end
end


max_peaks = max(all_0spikes46(3:158, :));
min_peaks = min(all_0spikes46(3:158, :));
include = find(max_peaks < 0.15 & min_peaks > -0.15)
avg_0spikes46 = mean(all_0spikes46(:,include),2);
std_0spikes46 = std(all_0spikes46(:,include)');
confint_top = avg_0spikes46 + 2* std_0spikes46'
confin_bot = avg_0spikes46 - 2 * std_0spikes46'

figure;
hold on
plot(all_0spikes46(3:158, include), 'Color', [0.5 0.5 0.5])
plot(avg_0spikes46(3:158), 'k', 'LineWidth', 3)
plot(confint_top(3:158), 'k', 'LineWidth', 3)
plot(confin_bot(3:158), 'k', 'LineWidth', 3)
plot([3 158], [0.1 0.1], 'm', 'LineWidth', 2)
hold off
xlim([3 158])
ax=gca;
    xsize = length(avg_0spikes46);
    ax.XTick = [0, xsize/7, (xsize/7)*2, (xsize/7)*3, (xsize/7)*4, (xsize/7)*5, (xsize/7)*6, (xsize/7)*7];
    ax.XTickLabel = {'0','1', '2', '3', '4', '5', '6', '7'};
    xlabel('time (sec)')
    ylabel('\DeltaF/F_{0}')
    set(gca, 'FontSize', 30)
%% Get all 1 spike data together

counter = 1
for i = 1:length(CellsAll46)
    if CellsAll46(i).spikeCount == 1
        all_1spikes46(:,counter) = CellsAll46(i).trace(1:159);
        counter = counter + 1
    end
end
avg_1spikes46 = mean(all_1spikes46,2);
avg_1spikes46 = mean(all_1spikes46(:,include1),2);

max_peaks1 = max(all_1spikes46(3:158, :));
min_peaks1 = min(all_1spikes46(3:158, :));

include1 = find(max_peaks1 < 0.15 & min_peaks > -0.15)
figure;
hold on
plot(all_1spikes46(3:158, :), 'Color', [0.5 0.5 0.5])
plot(avg_1spikes46(3:158), 'k', 'LineWidth', 3)

%% Find source of very large peaks in 0 spikes data
ids = [];
for i = 1:length(CellsAll49)
    if CellsAll49(i).spikeCount == 0
        if max(CellsAll49(i).trace(3:158)) > 0.1
            ids = [ids; max(CellsAll49(i).trace(3:158)), CellsAll49(i).cellnum, CellsAll49(i).trialnum];
        end
    end
end
hist(ids(:,1))

ids46 = [];
for i = 1:length(CellsAll46)
    if CellsAll46(i).spikeCount == 0
        if max(CellsAll46(i).trace(3:158)) > 0.1
            ids46 = [ids46; max(CellsAll46(i).trace(3:158)), CellsAll46(i).cellnum, CellsAll46(i).trialnum];
        end
    end
end
hist(ids46(:,1))












