%% Dissertation figure 6: high corr ROIs

%% A - example high corr
% use tad 36 (t=22) lag0R MS to show
% highcorr = 2 (actual ROI = 
% nothighcorr = 10 (actual ROI = 

plot_xcorr(allData{1,22}.lag0R_MS_sorted([2,10], :), 'hot')


%% C/D - how many high corr ROIs?
% to fill in venn diagram
highcorr_list = zeros(length(allData),9);
for t = 1:length(allData)
    highcorr_list(t,1) = length(allData{1,t}.resp_ROIs); %all data
    if isfield(allData{1,t}, 'uniqueHighCorrROI_MS')
        highcorr_list(t,3) = length(allData{1,t}.uniqueHighCorrROI_MS); %MS
        if isfield(allData{1,t}, 'uniqueHighCorrROI_V')
            highcorr_list(t,4) = length(allData{1,t}.uniqueHighCorrROI_V); %V
            if isfield(allData{1,t}, 'uniqueHighCorrROI_M')
                tmp = [allData{1,t}.uniqueHighCorrROI_MS; allData{1,t}.uniqueHighCorrROI_V; allData{1,t}.uniqueHighCorrROI_M]
                highcorr_list(t,2) = length(unique(tmp)); %all high corr ROIs
                highcorr_list(t,5) = length(allData{1,t}.uniqueHighCorrROI_M); %M
                
            end
        end
    end
    if isfield(allData{1,t}, 'highcorr_MS_V') 
        %if isfield(allData{1,t}, 'highcorr_V')
            highcorr_list(t,7) = length(allData{1,t}.highcorr_MS_V); %MS and V overlap  
    end
        if isfield(allData{1,t}, 'highcorr_M_MS')
            highcorr_list(t,8) = length(allData{1,t}.highcorr_M_MS); %MS and V overlap
        end
   % end
    if isfield(allData{1,t}, 'highcorr_M_V') %& isfield(allData{1,t}, 'highcorr_V')
        highcorr_list(t,6) = length(allData{1,t}.highcorr_M_V); %M and V overlap
    end
    if isfield(allData{1,t}, 'highcorr_M_V_MS')
        highcorr_list(t,9) = length(allData{1,t}.highcorr_M_V_MS); %MS, M amd V overlap
    end
end

% Summary data (Fig 6D)
highcorr_sums(1,:) = sum(highcorr_list(s46_tads,:))
highcorr_sums(2,:) = sum(highcorr_list(s49_tads,:))

%% REcalculate by tadpole - use to generate summary venn diagram
for t = 1:size(highcorr_list, 1)
    HC_prop(t,:) = highcorr_list(t, 2:end) ./ highcorr_list(t,1);
end

HC__prop_avg(1,:) = nanmean(HC_prop(st46,:), 1)
HC__prop_avg(2,:) = nanmean(HC_prop(st49,:), 1)


%% B: proportion highcorr by stage and modality

for h = 2:size(highcorr_list,1)
    high_corr_prop(h,:) = highcorr_list(h,2:end) / highcorr_list(h,1)
end

% Mean, and stats for each
for r = 1:size(high_corr_prop,2)
    stats_propHC(r,1) = nanmean(high_corr_prop(st49,r));
    stats_propHC(r,2) = nanmean(high_corr_prop(st46,r));
    stats_propHC(r,3) = nanstd(high_corr_prop(st49,r));
    stats_propHC(r,4) = nanstd(high_corr_prop(st46,r));    
end


% use st46_X (all 1) and st49_X (all 2) and s49_tads and s46_tads
st46_X = ones(1, length(s46_tads))
st49_X = 2*ones(1, length(s49_tads))
m=15
st = 1
figure;
hold on
plot(st46_X, high_corr_prop(s46_tads, 2), 'oy', 'MarkerSize', m)%46 MS 
plot(st46_X+0.3, high_corr_prop(s46_tads, 3), 'or', 'MarkerSize', m) %46 V
plot(st46_X+0.6, high_corr_prop(s46_tads, 4), 'ob', 'MarkerSize', m) %46 MS
plot(st49_X, high_corr_prop(s49_tads, 2), 'oy', 'MarkerSize', m) %49 MS 
plot(st49_X+0.3, high_corr_prop(s49_tads, 3), 'or', 'MarkerSize', m) %49 V
plot(st49_X+0.6, high_corr_prop(s49_tads, 4), 'ob', 'MarkerSize', m) %49 MS

plot([(st-0.1) (st+0.1)], [stats_propHC(2,2), stats_propHC(2,2)], 'k', 'LineWidth', 6) %46 MS
plot([(st-0.1+0.3) (st+0.1+0.3)], [stats_propHC(3,2), stats_propHC(3,2)], 'k', 'LineWidth', 6) % 46 V
plot([(st-0.1+0.6) (st+0.1+0.6)], [stats_propHC(4,2), stats_propHC(4,2)], 'k', 'LineWidth', 6) % 46 V

plot([(st-0.1+1) (st+0.1+1)], [stats_propHC(2,1), stats_propHC(2,1)], 'k', 'LineWidth', 6) %46 MS
plot([(st-0.1+1.3) (st+0.1+1.3)], [stats_propHC(3,1), stats_propHC(3,1)], 'k', 'LineWidth', 6) % 46 V
plot([(st-0.1+1.6) (st+0.1+1.6)], [stats_propHC(4,1), stats_propHC(4,1)], 'k', 'LineWidth', 6) % 46 V

hold off 
xlim([0.6 3])
ax = gca;
ax.XTick = [1 1.3 1.6 2 2.3 2.6]
ax.XTickLabel = {'MS', 'V', 'M', 'MS', 'V', 'M'}
ylabel('proportion of ROIs')
set(gca, 'FontSize', 30)
saveas(gcf, 'prop ROIs highcorr by mod and st', 'png')
saveas(gcf, 'prop ROIs highcorr by mod and st', 'epsc2')

% Stats
high_corr_prop_stats(1,1) = nanmean(high_corr_prop(s46_tads, 2))
high_corr_prop_stats(2,1) = nanmean(high_corr_prop(s46_tads, 3))
high_corr_prop_stats(3,1) = nanmean(high_corr_prop(s46_tads, 4))
high_corr_prop_stats(4,1) = nanmean(high_corr_prop(s49_tads, 2))
high_corr_prop_stats(5,1) = nanmean(high_corr_prop(s49_tads, 3))
high_corr_prop_stats(6,1) = nanmean(high_corr_prop(s49_tads, 4))

high_corr_prop_stats(1,2) = nanstd(high_corr_prop(s46_tads, 2))
high_corr_prop_stats(2,2) = nanstd(high_corr_prop(s46_tads, 3))
high_corr_prop_stats(3,2) = nanstd(high_corr_prop(s46_tads, 4))
high_corr_prop_stats(4,2) = nanstd(high_corr_prop(s49_tads, 2))
high_corr_prop_stats(5,2) = nanstd(high_corr_prop(s49_tads, 3))
high_corr_prop_stats(6,2) = nanstd(high_corr_prop(s49_tads, 4))


%% E (Figure 7): any differences between highcorr and nothighcorr?

% Use allRespROIs to assemble info
% add a field to allData{1,t} that is all unique highcorr ROIs
for t = 1:length(allData)
    tmp = [];
    if isfield(allData{1,t}, 'uniqueHighCorrROI_MS') & isfield(allData{1,t}, 'uniqueHighCorrROI_V') & isfield(allData{1,t}, 'uniqueHighCorrROI_M')
        allData{1,t}.allhighcorrROI = unique([allData{1,t}.uniqueHighCorrROI_MS; allData{1,t}.uniqueHighCorrROI_V; allData{1,t}.uniqueHighCorrROI_M]);
    else 
        if isfield(allData{1,t}, 'uniqueHighCorrROI_MS')
            tmp = [tmp; allData{1,t}.uniqueHighCorrROI_MS];
        end
        if isfield(allData{1,t}, 'uniqueHighCorrROI_V')
            tmp = [tmp; allData{1,t}.uniqueHighCorrROI_V];
        end
        if isfield(allData{1,t}, 'uniqueHighCorrROI_M')
            tmp = [tmp; allData{1,t}.uniqueHighCorrROI_M];
        end
        tmp1 = unique(tmp)
        allData{1,t}.allhighcorrROI = tmp1;
    end
end
          
% add a column to allRespROIs that defines highcorr = 1  and nothighcorr = 0
for r = 1:size(allRespROIs,1)
    if ismember(allRespROIs(r,2), allData{1,allRespROIs(r,1)}.allhighcorrROI)
        allRespROIs(r, 30) = 1;
    else
        allRespROIs(r, 30) = 0;
    end
end

st46r = find(allRespROIs(:,3) == 46)
st49r = find(allRespROIs(:,3) == 49)
% index by stage (col 28) and hc/nhc (col 30) to group and stat test
hc = find(allRespROIs(:,30) == 1);
nhc = find(allRespROIs(:,30) == 0);
st46h = intersect(st46r, hc)
st49h = intersect(st49r, hc)
st46n = intersect(st46r, nhc)
st49n = intersect(st49r, nhc)

%labels = {'area MS'; 'area V'; 'area M'; 'area NS'; 'peak MS'; 'peak V'; 'peak M'; 'peak NS'; ... 
%    'MSIndex peak'; 'MSIndex area'; 'MSIndex onsettime'; 'nothing'; 'MSEnh peak'; 'unimax peak'; 'unimax stimtype'; ... 
%    'onset time MS'; 'onset time V'; 'onset time M'; 'onset time NS';...
%    'onset time SD MS'; 'onset time SD V'; 'onset time SD M'; 'onset time SD NS'; ...
%     'MSEnh num resp'; 'Uni bias num resp'};

% Stat test hc vs nhc for each stage seperately (b/c no 2 way ANOVA for
% nonparametric data)
for i = 3:27
    hc = allRespROIs(st46h, i);
    nhc = allRespROIs(st46h, i);
    allhighcorr46_H(i-2) = kstest2(hc, nhc);
end
diff_vars = labels(find(allhighcorr46_H))
% no stat diffs
for i = 3:27
    hc = allRespROIs(st49h, i);
    nhc = allRespROIs(st49h, i);
    allhighcorr49_H(i-2) = kstest2(hc, nhc);
end
diff_vars = labels(find(allhighcorr49_H))
% no stat diffs

for i = 3:27
    [p{i},tbl{i},stats{i}] = anovan(allRespROIs(:,i), {allRespROIs(:,28), allRespROIs(:,30)}, 'varnames', {'stage', 'HC status'})
end
% confusing

% eliminate stage as variable
for i = 3:27
    hc_vals = allRespROIs(hc, i);
    nhc_vals = allRespROIs(nhc, i);
    [allhighcorr_H(i-2, 1), allhighcorr_H(i-2, 2), allhighcorr_H(i-2, 3)] = kstest2(hc_vals, nhc_vals);
end
diff_vars = labels(find(allhighcorr_H(:,1)))

%% Make ECDF plots for variables that are significant

% used this image to select colors: http://www.somersault1824.com/wp-content/uploads/2015/02/color-blindness-palette.png

for i = 1:size(allhighcorr_H, 1)
    if allhighcorr_H(i, 1)
        [f1, x1] = ecdf(allRespROIs(hc, i+2)); 
        [f2, x2] = ecdf(allRespROIs(nhc, i+2)); 
        figure;
        hold on
        plot(x1, f1, 'Color', [0 .29, .29], 'LineWidth', 3) % hc in blue
        plot(x2, f2, 'Color', [.71 .29, 1], 'LineWidth', 3) %nhc in purple
        hold off
        ylabel('ROIs')
        xlabel(labels{i})
        %xlim([-6 15])
        set(gca,'FontSize',20)
        fig_filename = sprintf('ecdf hc_nhc %s', labels{i})
        saveas(gcf, fig_filename, 'png')
        saveas(gcf, fig_filename, 'epsc2')
        close;
    end
end

% make histograms

for i = 1:size(allhighcorr_H, 1)
    if allhighcorr_H(i, 1)
        figure;
        histogram(allRespROIs(hc, i+2), 40, 'FaceColor', [.71 .29, 1], 'EdgeColor', [.71 .29, 1]) %nhc in purple

        ylabel('ROI Count')
        xlabel(labels{i})
        %xlim([-6 15])
        set(gca,'FontSize',20)
        fig_filename = sprintf('hist hc_nhc NHC %s', labels{i})
        saveas(gcf, fig_filename, 'png')
        saveas(gcf, fig_filename, 'epsc2')
        close;
    end
end

for i = 1:size(allhighcorr_H, 1)
    if allhighcorr_H(i, 1)
        figure;

        histogram(allRespROIs(hc, i+2), 40, 'FaceColor', [0 .29, .29], 'EdgeColor', [0 .29, .29]) % hc in blue

        ylabel('ROI Count')
        xlabel(labels{i})
        %xlim([-6 15])
        set(gca,'FontSize',20)
        fig_filename = sprintf('hist hc_nhc HC %s', labels{i})
        saveas(gcf, fig_filename, 'png')
        saveas(gcf, fig_filename, 'epsc2')
        close;
    end
end




%% Make specific histograms for figure:

Labels = {'onset time (sec)', 'SD of onset time (sec)', 'peak (\DeltaF/F_{0})'}
filenames = {'onset time MS USE', 'onset time SD MS USE', 'unimax peak USE'}
list = [14, 16, 20];
limits = [-0.05, 2; 0, 6; 0, 3]
for i = 1:length(list)
        figure;
        histogram(allRespROIs(hc, list(i)+2), 40, 'FaceColor', [0 .29, .29], 'EdgeColor', [0 .29, .29]) % hc in blue
        ylabel('ROI Count')
        xlabel(Labels{i})
        xlim(limits(i, :))
        set(gca,'FontSize',20)
        fig_filename = sprintf('hist hc_nhc HC %s', filenames{i})
        saveas(gcf, fig_filename, 'png')
        saveas(gcf, fig_filename, 'epsc2')
        %close;
        figure;
        histogram(allRespROIs(nhc, list(i)+2), 40, 'FaceColor', [.71 .29, 1], 'EdgeColor', [.71 .29, 1]) %nhc in purple
        ylabel('ROI Count')
        xlabel(Labels{i})
        xlim(limits(i, :))
        set(gca,'FontSize',20)
        fig_filename = sprintf('hist hc_nhc NHC %s', filenames{i})
        saveas(gcf, fig_filename, 'png')
        saveas(gcf, fig_filename, 'epsc2')
        %close;
        
end

for i = 1:length(list)
        [f1, x1] = ecdf(allRespROIs(hc, list(i)+2)); 
        [f2, x2] = ecdf(allRespROIs(nhc, list(i)+2)); 
        figure;
        hold on
        plot(x1, f1, 'Color', [0 .29, .29], 'LineWidth', 3) % hc in blue
        plot(x2, f2, 'Color', [.71 .29, 1], 'LineWidth', 3) %nhc in purple
        hold off
        ylabel('ROI proportion')
        xlabel(Labels{i})
        xlim(limits(i, :))
        set(gca,'FontSize',20)
        fig_filename = sprintf('ecdf hc_nhc NHC %s', filenames{i})
        saveas(gcf, fig_filename, 'png')
        saveas(gcf, fig_filename, 'epsc2')
end

%% Notes
% make ECDF plots for each var 
for i = 3:27
    s46 = allRespROIs(st46, i);
    s49 = allRespROIs(st49, i);
    figure;
    hold on
    ecdf(s46)

    %set(h, 'Color', 'g')
    ecdf(s49)
    h = get(gca, 'children')
    set(h, 'LineWidth', 3)
    %j = get(gca, 'children')
    %set(j, 'LineWidth', 3)
    %set(j, 'Color', [0.5 0 0.5])
    hold off
    title(sprintf('%s ECDF by stage', labels{i-2}))
    xlabel(labels{i-2})
    ylabel('ROI count')
    fig_filename = sprintf('st 46 vs 49 ECDF of %s', labels{i-2})
    saveas(gcf, fig_filename, 'png')
    close;
end

% are any statistically different? (using kstest2)
for i = 3:27
    s46 = allRespROIs(st46, i);
    s49 = allRespROIs(st49, i);
    allRespROIs_H(i-2) = kstest2(s46, s49);
end
diff_vars = labels(find(allRespROIs_H))


