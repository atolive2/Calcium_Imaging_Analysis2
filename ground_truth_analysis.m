% cell attached analysis

%Plot df_f0 of patched cell
hold on
for j = 1:size(tadpole.df_f0,1)
    plot(tadpole.df_f0{1,j})
end
hold off
ax=gca;
xsize = length(tadpole.df_f0{1,j});
ax.XTick = [0, xsize/7, (xsize/7)*2, (xsize/7)*3, (xsize/7)*4, (xsize/7)*5, (xsize/7)*6, (xsize/7)*7];
ax.XTickLabel = {'0','1', '2', '3', '4', '5', '6', '7'};
title(sprintf('Exp %d cell 1 all trials', tadpole.expnum));
xlabel('time(s)');
ylabel('\Delta F/F_{0}');


% smooth df_f0
for i = 1:size(tadpole.df_f0, 1)
    for j = 1:size(tadpole.df_f0, 2)
        tadpole.filtered{i,j} = smooth(tadpole.df_f0{i,j}(:,:), 8, 'moving');
    end
end
% plot smoothed data
figure;
hold on
for j = 1:size(tadpole.filtered,1)
    plot(tadpole.filtered{1,j})
end
hold off
ax=gca;
xsize = length(tadpole.filtered{1,j});
ax.XTick = [0, xsize/7, (xsize/7)*2, (xsize/7)*3, (xsize/7)*4, (xsize/7)*5, (xsize/7)*6, (xsize/7)*7];
ax.XTickLabel = {'0','1', '2', '3', '4', '5', '6', '7'};
title(sprintf('Exp %d cell 1 all trials, smoothed', tadpole.expnum));
xlabel('time(s)');
ylabel('\Delta F/F_{0}');
