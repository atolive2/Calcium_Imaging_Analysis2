%% Plot peak of exp 5 and then all exps

% assuming you have opened exp5tadpole

%subplot 1: peak - lgst uni vs MSenh all
% get values for plotting Y=X line
% uni_mm = [min(min([peak_V(1,:) peak_M(1,:)])), max(max([peak_V(1,:), peak_M(1,:)]))];
% multi_mm = [min(min([peak_V(2,:), peak_M(2,:)])), max(max([peak_V(2,:), peak_M(2,:)]))];
% Y = (min([uni_mm, multi_mm])):0.01:(max(([uni_mm, multi_mm])));

subplot(1,2,1)
hold on
plot(peak_V(1,:), peak_V(2,:), 'ro', 'MarkerFaceColor', 'r')
plot(peak_M(1,:), peak_M(2,:), 'bo', 'MarkerFaceColor', 'b')
% plot(Y, Y, 'b-')
ylim([0 0.6])
xlim([0 0.6])
title('Largest Unisensory Response vs Multisensory Enhancement')
xlabel('largest unisensory mean')
ylabel('multisensory enhancement mean')
legend('Vis', 'Mech', 'location', 'best')
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
ylim([0 0.6])
xlim([0 0.6])
title('Linear Sum vs Multisensory Response')
xlabel('linear sum mean')
ylabel('multisensory response mean')
legend('Vis', 'Mech', 'Y=X', 'location', 'best')
hold off

%% Find ROIs that demostrate PIE

% exp 7 (struct9), roi 10. 
figure;
    for i = 1:size(tadpole{1,9}.df_f0,2)
        hold on
        if tadpole{1,9}.stimorder(1, i) == 1
            plot(tadpole{1,9}.df_f0{10,i}, 'color', [0.5 0 0.5])
        elseif tadpole{1,9}.stimorder(1, i) == 2
            plot(tadpole{1,9}.df_f0{10,i}, 'b')
        elseif tadpole{1,9}.stimorder(1, i) == 3
            plot(tadpole{1,9}.df_f0{10,i}, 'r')
        elseif tadpole{1,9}.stimorder(1, i) == 4
            plot(tadpole{1,9}.df_f0{10,i}, 'color', [0.5 0.5 0.5])
        end
    end
% hahaha nope. 

% Smooth all ROIs over all exps, and then plot all like above and look for
% good ones.
% apply filter to all trials (moving-average filter)
for t = 1:length(tadpole)
    for i = 1:size(tadpole{1,t}.df_f0, 1)
        for j = 1:size(tadpole{1,t}.df_f0, 2)
            tadpole{1,t}.filtered{i,j} = smooth(tadpole{1,t}.df_f0{i,j}(:,:), 8, 'moving');
        end
    end
end

% plot smoothed for all ROIs, all tads
for t = 1:length(tadpole) % over all tads
    for r = 1:size(tadpole{1,t}.filtered,1) % over all rois
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
        title(sprintf('Exp %d ROI %d smoothed all', tadpole{1,t}.expnum, r))
        xlabel('time (s)')
        ylabel('\DeltaF/F_{0}')
        ax=gca;
        xsize = length(tadpole{1,t}.filtered{1,1});
        ax.XTick = [0, xsize/7, (xsize/7)*2, (xsize/7)*3, (xsize/7)*4, (xsize/7)*5, (xsize/7)*6, (xsize/7)*7];
        ax.XTickLabel = {'0','1', '2', '3', '4', '5', '6', '7'};
        hold off
          fig_filename=sprintf('F:/Calcium_Imaging_Analysis/tadpoles_byexp/figures/smoothedROIs/exp%droi%dsmoothed_all.png', tadpole{1,t}.expnum, r)
        saveas(gcf,fig_filename,'png');
        close;
        clear('fig_filename')
    end
end












