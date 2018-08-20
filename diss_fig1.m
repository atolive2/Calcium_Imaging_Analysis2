%% For methods examples
% ROI 3
% plot all MS trials on top of each other
MSstims = [1 5 8 9]
Vstims = [2 6]
Mstims = [3 7]
t = 20 %tad 34 is in allData{1,20}
roi = 3
figure;
hold on
for s = 1:length(allData{1,t}.stimorder)
    if ismember(allData{1,t}.stimorder(s), MSstims)
        if max(allData{1,t}.df_f0{roi,s}) < 1
            if min(allData{1,t}.df_f0{roi,s}) > -0.1
                plot(allData{1,t}.df_f0{roi,s}, 'LineWidth', 1)
            end
        end
    end
end
hold off
    ax=gca;
    xsize = length(allData{1,t}.df_f0{roi,s});
    ax.XTick = [0, xsize/7, (xsize/7)*2, (xsize/7)*3, (xsize/7)*4, (xsize/7)*5, (xsize/7)*6, (xsize/7)*7];
    ax.XTickLabel = {'0','1', '2', '3', '4', '5', '6', '7'};
    ylabel('\DeltaF/F_{0}', 'fontsize', 30)
    xlabel('time (sec)', 'fontsize', 30)
    set(gca, 'fontsize', 20)
    fig_filename = sprintf('Tad 34 ROI %d all good MS trials df_f0)', roi)
    saveas(gcf, fig_filename, 'epsc2')
    
% ROI 63
roi = 63
figure;
hold on
for s = 1:length(allData{1,t}.stimorder)
    if ismember(allData{1,t}.stimorder(s), MSstims)
        if max(allData{1,t}.df_f0{roi,s}) < 1
            if min(allData{1,t}.df_f0{roi,s}) > -0.1
                plot(allData{1,t}.df_f0{roi,s}, 'LineWidth', 1)
            end
        end
    end
end
hold off
    ax=gca;
    xsize = length(allData{1,t}.df_f0{roi,s});
    ax.XTick = [0, xsize/7, (xsize/7)*2, (xsize/7)*3, (xsize/7)*4, (xsize/7)*5, (xsize/7)*6, (xsize/7)*7];
    ax.XTickLabel = {'0','1', '2', '3', '4', '5', '6', '7'};
    ylabel('\DeltaF/F_{0}', 'fontsize', 30)
    xlabel('time (sec)', 'fontsize', 30)
    set(gca, 'fontsize', 20)
    fig_filename = sprintf('Tad 34 ROI %d all good MS trials df_f0)', roi)
    saveas(gcf, fig_filename, 'epsc2')