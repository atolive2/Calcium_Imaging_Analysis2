%% analysis of Stage 46 tadpoles

% Exp 32 only
% each trial, all ROIs
for i = 1:size(tadpole.df_f0,2)
    figure;
    hold on
    for j = 1:size(tadpole.df_f0,1)
        plot(tadpole.df_f0{j,i})
    end
    ax=gca;
    xsize = length(tadpole.df_f0{j,i});
    ax.XTick = [0, xsize/7, (xsize/7)*2, (xsize/7)*3, (xsize/7)*4, (xsize/7)*5, (xsize/7)*6, (xsize/7)*7];
    ax.XTickLabel = {'0','1', '2', '3', '4', '5', '6', '7'};
    title(sprintf('Exp %d trial %d all ROIs', tadpole.expnum, i));
    xlabel('time(s)');
    ylabel('\DeltaF/F_{0}');
    hold off
    fig_filename=sprintf([tadpole.figure_filepath 'exp%dtrial%d_dff0.png'], tadpole.expnum, i);
    saveas(gcf,fig_filename,'png');
    close;
    clear('fig_filename')
end

% Each ROI, all trials

%%%%%%% need to eliminate bad trials efficiently%%%%%%%%%%%%%%%%%%%%%%%%%%%
trials = 1:size(tadpole.df_f0,2) %standard - keep all trials
if isfield(tadpole, 'badtrials') % if there are bad trials
    for a = 1:length(tadpole.badtrials)
        trials(trials == tadpole.badtrials(a)) = []
    end
    
end
length(trials) %should be size(tadpole.df_f0) - length(tadpole.badtrials) e.g. the number of good trials

for i = 1:size(tadpole.df_f0,1)
    figure;
    hold on
    for j = 1:length(trials)
        plot(tadpole.df_f0{i,trials(j)})
    end
    ax=gca;
    xsize = length(tadpole.df_f0{i,j});
    ax.XTick = [0, xsize/7, (xsize/7)*2, (xsize/7)*3, (xsize/7)*4, (xsize/7)*5, (xsize/7)*6, (xsize/7)*7];
    ax.XTickLabel = {'0','1', '2', '3', '4', '5', '6', '7'};
    title(sprintf('Exp %d ROI %d all Trials', tadpole.expnum, i));
    xlabel('time(s)');
    ylabel('\DeltaF/F_{0}');
    hold off
    fig_filename=sprintf([tadpole.figure_filepath 'exp%dROI%d_dff0.png'], tadpole.expnum, i);
    saveas(gcf,fig_filename,'png');
    close;
    clear('fig_filename')
end