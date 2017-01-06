%% Using exp 5, generate ROI info for PIE figure 1 (all traces for an ROI)
% assumes you have opened the exp 5 file

% find the right ROIs
% ROI 5 - multi
% ROI 119 - vis
% ROI 63 - mech

% plot all trials for that ROI, with multi = purple, vis = red, mech =
% blue, no stim = black

% ROI 5
figure;
for i = 1:size(tadpole.df_f0, 2)
    hold on
    if tadpole.stimorder(1, i) == 1
        plot(tadpole.df_f0{5,i}, 'color', [0.5 0 0.5])
    elseif tadpole.stimorder(1, i) == 2
        plot(tadpole.df_f0{5,i}, 'b')
    elseif tadpole.stimorder(1, i) == 3
        plot(tadpole.df_f0{5,i}, 'r')
    elseif tadpole.stimorder(1, i) == 4
        plot(tadpole.df_f0{5,i}, 'color', [0.5 0.5 0.5])
    end
end
title('Exp 5 ROI 5')
xlabel('time (s)')
ylabel('\DeltaF/F_{0}')
ax=gca;
xsize = length(tadpole.df_f0{5,1});
ax.XTick = [0, xsize/7, (xsize/7)*2, (xsize/7)*3, (xsize/7)*4, (xsize/7)*5, (xsize/7)*6, (xsize/7)*7];
ax.XTickLabel = {'0','1', '2', '3', '4', '5', '6', '7'};

% ROI 119
figure;
for i = 1:size(tadpole.df_f0, 2)
    hold on
    if tadpole.stimorder(1, i) == 1
        plot(tadpole.df_f0{117,i}, 'color', [0.5 0 0.5])
    elseif tadpole.stimorder(1, i) == 2
        plot(tadpole.df_f0{117,i}, 'b')
    elseif tadpole.stimorder(1, i) == 3
        plot(tadpole.df_f0{117,i}, 'r')
    elseif tadpole.stimorder(1, i) == 4
        plot(tadpole.df_f0{117,i}, 'color', [0.5 0.5 0.5])
    end
end
title('Exp 5 ROI 119')
xlabel('time (s)')
ylabel('\DeltaF/F_{0}')
ax=gca;
xsize = length(tadpole.df_f0{5,1});
ax.XTick = [0, xsize/7, (xsize/7)*2, (xsize/7)*3, (xsize/7)*4, (xsize/7)*5, (xsize/7)*6, (xsize/7)*7];
ax.XTickLabel = {'0','1', '2', '3', '4', '5', '6', '7'};

% ROI 138
figure;
for i = 1:size(tadpole.df_f0, 2)
    hold on
    if tadpole.stimorder(1, i) == 1
        plot(tadpole.df_f0{62,i}, 'color', [0.5 0 0.5])
    elseif tadpole.stimorder(1, i) == 2
        plot(tadpole.df_f0{62,i}, 'b')
    elseif tadpole.stimorder(1, i) == 3
        plot(tadpole.df_f0{62,i}, 'r')
    elseif tadpole.stimorder(1, i) == 4
        plot(tadpole.df_f0{62,i}, 'color', [0.5 0.5 0.5])
    end
end
title('Exp 5 ROI 63')
xlabel('time (s)')
ylabel('\DeltaF/F_{0}')
ax=gca;
xsize = length(tadpole.df_f0{5,1});
ax.XTick = [0, xsize/7, (xsize/7)*2, (xsize/7)*3, (xsize/7)*4, (xsize/7)*5, (xsize/7)*6, (xsize/7)*7];
ax.XTickLabel = {'0','1', '2', '3', '4', '5', '6', '7'};