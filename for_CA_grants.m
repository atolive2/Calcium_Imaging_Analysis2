%% figures for CA grant (10/24/2016)
% gcamp ours

for i = 1:size(tadpole.df_f0,1)
    figure;
    hold on
    for j = 1:size(tadpole.df_f0,2)
        plot(tadpole.df_f0{i,j})
    end
    title(sprintf('ROI %d \DeltaF/F_{0}', i))
    ax=gca;
    xsize = length(tadpole.df_f0{i,j});
    ax.XTick = [0, xsize/7, (xsize/7)*2, (xsize/7)*3, (xsize/7)*4, (xsize/7)*5, (xsize/7)*6, (xsize/7)*7];
    ax.XTickLabel = {'0','1', '2', '3', '4', '5', '6', '7'};
    xlabel('time(s)');
    ylabel('\DeltaF/F_{0}');
    hold off
    fig_filename=sprintf([tadpole.figure_filepath 'ROI%d_df_F0.png'],i);
    saveas(gcf,fig_filename,'png');
    close;
    clear('fig_filename')
end

hold on
for i = 1:size(tadpole.df_f0,2)
    figure;
    hold on
    for j = 1:(size(tadpole.df_f0,1)-1)
        plot(tadpole.df_f0{j,i})
    end
    ax=gca;
    xsize = length(tadpole.df_f0{j,i});
    ax.XTick = [0, xsize/7, (xsize/7)*2, (xsize/7)*3, (xsize/7)*4, (xsize/7)*5, (xsize/7)*6, (xsize/7)*7];
    ax.XTickLabel = {'0','1', '2', '3', '4', '5', '6', '7'};
    title(sprintf('Exp %d trial %d all ROIs \DeltaF/F_{0}', tadpole.expnum, i));
    xlabel('time(s)');
    ylabel('raw pixel intensity');
    hold off
    fig_filename=sprintf([tadpole.figure_filepath 'exp%dtrial%d_df_f0.png'], tadpole.expnum, i);
    saveas(gcf,fig_filename,'png');
    close;
    clear('fig_filename')
end

% plot just trial 1, rois 2,3,4,5,6
roirng = [2 3 4 5 6];
figure;
    hold on
    for j = 1:length(roirng)
        plot(tadpole.df_f0{roirng(j),i})
    end
    ax=gca;
    xsize = length(tadpole.df_f0{j,i});
    ax.XTick = [0, xsize/7, (xsize/7)*2, (xsize/7)*3, (xsize/7)*4, (xsize/7)*5, (xsize/7)*6, (xsize/7)*7];
    ax.XTickLabel = {'0','1', '2', '3', '4', '5', '6', '7'};
    title(sprintf('Exp %d trial %d all ROIs \DeltaF/F_{0}', tadpole.expnum, i));
    xlabel('time(s)');
    ylabel('\DeltaF/F_{0}');
    legend('ROI 1', 'ROI 2', 'ROI 3', 'ROI 4', 'ROI 5')
    hold off
    fig_filename=sprintf([tadpole.figure_filepath 'exp%dtrial%d_df_f0.png'], tadpole.expnum, i);
    saveas(gcf,fig_filename,'png');
    close;
    clear('fig_filename')
    
    
%% OGB1 Examples

% apply filter to all trials (moving-average filter)
for i = 1:size(tadpole.df_f0, 1)
    for j = 1:size(tadpole.df_f0, 2)
        filtered{i,j} = smooth(tadpole.df_f0{i,j}(:,:), 8, 'moving');
    end
end

% plot blocks 1 and 2 by ROIs of interest
rois = [18 19 31 32 44 54];


for r = 1:length(rois)
    figure;
    for i = 1:size(filtered(rois(r),1:24), 2)
        hold on
        if tadpole.stimorder(1, i) == 1
            plot(filtered{rois(r),i}, 'color', [0.5 0 0.5])
        elseif tadpole.stimorder(1, i) == 2
            plot(filtered{rois(r),i}, 'b')
        elseif tadpole.stimorder(1, i) == 3
            plot(filtered{rois(r),i}, 'r')
        elseif tadpole.stimorder(1, i) == 4
            plot(filtered{rois(r),i}, 'color', [0.5 0.5 0.5])
        end
    end
    title(sprintf('Exp 6 ROI %d', rois(r)))
    xlabel('time (s)')
    ylabel('\DeltaF/F_{0}')
    ax=gca;
    xsize = 220
    %length(tadpole.df_f0{,1});
    ax.XTick = [0, xsize/7, (xsize/7)*2, (xsize/7)*3, (xsize/7)*4, (xsize/7)*5, (xsize/7)*6, (xsize/7)*7];
    ax.XTickLabel = {'0','1', '2', '3', '4', '5', '6', '7'};
    hold off
%     fig_filename=sprintf([tadpole.figure_filepath 'exp%dtrial%d.png'], tadpole.expnum, i);
%     saveas(gcf,fig_filename,'png');
%    close;
%    clear('fig_filename')
end

for r = 1:length(rois)
    for i = 1:24
    peak(r,i) = max(filtered{r,i})
    end
end

% find exactly the right traces and plot. 
% all are 1 2 3 4 repeating order. (1,5,8,13,17,20 = 1)
multi = [1,5,8,13,17,20]
vis = [2,6,9,14,18,21]
mech=[3,7,10,15,19,23]
% what traces do you want?
figure;
hold on
plot(filtered{54,multi(1)}, 'y')
plot(filtered{54,multi(2)}, 'm')
plot(filtered{54,multi(3)}, 'c')
plot(filtered{54,multi(4)}, 'r')
plot(filtered{54,multi(5)}, 'g')
plot(filtered{54,multi(6)}, 'b')
title('multi')
hold off
figure;
hold on
plot(filtered{54,vis(1)}, 'y')
plot(filtered{54,vis(2)}, 'm')
plot(filtered{54,vis(3)}, 'c')
plot(filtered{54,vis(4)}, 'r')
plot(filtered{54,vis(5)}, 'g')
plot(filtered{54,vis(6)}, 'b')
title('vis')
hold off
figure;
hold on
plot(filtered{54,mech(1)}, 'y')
plot(filtered{54,mech(2)}, 'm')
plot(filtered{54,mech(3)}, 'c')
plot(filtered{54,mech(4)}, 'r')
plot(filtered{54,mech(5)}, 'g')
plot(filtered{54,mech(6)}, 'b')
title('mech')
hold off


% Exp 6 ROI 31
multi31 = [13,17,20]
vis31 = [6,9,18]
mech31 =[7,10,23]
%multi
figure;
hold on
for i = 1:length(multi31)
    plot(filtered{31,multi31(i)}, 'color', [0.5 0 0.5])
end
%vis
figure;
hold on
for i = 1:length(vis31)
    plot(filtered{31,vis31(i)}, 'b')
end
% mech
figure;
hold on
for i = 1:length(mech31)
    plot(filtered{31,mech31(i)}, 'r')
end
axis([2 219 -0.08 0.3])

% Exp 6 ROI 54
multi54 = [5,8,20]
vis54 = [6,18,21]
mech54 =[3,15,23]
%multi
figure;
hold on
for i = 1:length(multi54)
    plot(filtered{54,multi54(i)}, 'color', [0.5 0 0.5])
end
%vis
figure;
hold on
for i = 1:length(vis54)
    plot(filtered{54,vis54(i)}, 'b')
end
% mech
figure;
hold on
for i = 1:length(mech54)
    plot(filtered{54,mech54(i)}, 'r')
end
axis([2 219 -0.2 0.8])

% To get a primarily V ROI, move to exp 8, ROI 34. 
for i = 1:size(tadpole.df_f0, 1)
    for j = 1:size(tadpole.df_f0, 2)
        filtered{i,j} = smooth(tadpole.df_f0{i,j}(:,:), 8, 'moving');
    end
end

rois = 34
for r = 1:length(rois)
    figure;
    for i = 1:size(filtered(rois(r),1:12), 2)
        hold on
        if tadpole.stimorder(1, i) == 1
            plot(filtered{rois(r),i}, 'color', [0.5 0 0.5])
        elseif tadpole.stimorder(1, i) == 2
            plot(filtered{rois(r),i}, 'b')
        elseif tadpole.stimorder(1, i) == 3
            plot(filtered{rois(r),i}, 'r')
        elseif tadpole.stimorder(1, i) == 4
            plot(filtered{rois(r),i}, 'color', [0.5 0.5 0.5])
        end
    end
    title(sprintf('Exp 6 ROI %d', rois(r)))
    xlabel('time (s)')
    ylabel('\DeltaF/F_{0}')
    ax=gca;
    xsize = 220
    ax.XTick = [0, xsize/7, (xsize/7)*2, (xsize/7)*3, (xsize/7)*4, (xsize/7)*5, (xsize/7)*6, (xsize/7)*7];
    ax.XTickLabel = {'0','1', '2', '3', '4', '5', '6', '7'};
    hold off

end

multi = [1,5,8,13,17,21]
vis = [2,6,9,14,18,22]
mech=[3,7,10,15,19,23]
% what traces do you want?
figure;
hold on
plot(filtered{34,multi(1)}, 'y')
plot(filtered{34,multi(2)}, 'm')
plot(filtered{34,multi(3)}, 'c')
plot(filtered{34,multi(4)}, 'r')
plot(filtered{34,multi(5)}, 'g')
plot(filtered{34,multi(6)}, 'b')
title('multi')
hold off
figure;
hold on
plot(filtered{34,vis(1)}, 'y')
plot(filtered{34,vis(2)}, 'm')
plot(filtered{34,vis(3)}, 'c')
plot(filtered{34,vis(4)}, 'r')
plot(filtered{34,vis(5)}, 'g')
plot(filtered{34,vis(6)}, 'b')
title('vis')
hold off
figure;
hold on
plot(filtered{34,mech(1)}, 'y')
plot(filtered{34,mech(2)}, 'm')
plot(filtered{34,mech(3)}, 'c')
plot(filtered{34,mech(4)}, 'r')
plot(filtered{34,mech(5)}, 'g')
plot(filtered{34,mech(6)}, 'b')
title('mech')
hold off
axis([0 160 -0.05 0.2])

% Exp 8 ROI 34
multi34 = [1,5,8]
vis34 = [9,14,21]
mech34 = [3,7,23]
%multi
figure;
hold on
for i = 1:length(multi34)
    plot(filtered{34,multi34(i)}, 'color', [0.5 0 0.5])
end
%vis
figure;
hold on
for i = 1:length(vis34)
    plot(filtered{34,vis34(i)}, 'b')
end
% mech
figure;
hold on
for i = 1:length(mech34)
    plot(filtered{34,mech34(i)}, 'r')
end
axis([10 136 -0.05 0.2])
















