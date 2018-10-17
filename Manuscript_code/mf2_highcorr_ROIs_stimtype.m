%% Using analysis from manuscript_figs_2, assess correlations some other way

%% Are cells more correlated in response to MS or US? 

% plot MS maxR vs each US maxR

x_line = [0 1];
y_line = [0 1];
for k = 1:length(allData)
    if ~isempty(allData(k).respROIdff0_maxR_sq)
        figure;
        %MS vs V
        subplot(2,2,1)
        hold on
        for r = 1:size(allData(k).respROIdff0_maxR_sq, 2)
            plot(squeeze(allData(k).respROIdff0_maxR_sq(1,r,:)), squeeze(allData(k).respROIdff0_maxR_sq(2,r,:)), 'o');
        end
        plot(x_line, y_line, 'k', 'LineWidth', 2)
        axis tight
        title('MS vs V')
        xlabel('max R MS')
        ylabel('max R V')
        hold off 
        % MS vs M
        subplot(2,2,2)
        hold on
        for r = 1:size(allData(k).respROIdff0_maxR_sq, 2)
            plot(squeeze(allData(k).respROIdff0_maxR_sq(1,r,:)), squeeze(allData(k).respROIdff0_maxR_sq(3,r,:)), 'o');
        end
        plot(x_line, y_line, 'k', 'LineWidth', 2)
        axis tight
        title('MS vs M')
        xlabel('max R MS')
        ylabel('max R M')
        hold off 
        % MS vs N
        subplot(2,2,3)
        hold on
        for r = 1:size(allData(k).respROIdff0_maxR_sq, 2)
            plot(squeeze(allData(k).respROIdff0_maxR_sq(1,r,:)), squeeze(allData(k).respROIdff0_maxR_sq(4,r,:)), 'o');
        end
        plot(x_line, y_line, 'k', 'LineWidth', 2)
        axis tight
        title('MS vs N')
        xlabel('max R MS')
        ylabel('max R N')
        hold off 
        % overlaid histogram a la diss figure ??
        subplot(2,2,4) 
        hold on 
        histogram(reshape(allData(k).respROIdff0_maxR_sq(3,:,:), 1, []), 'FaceColor', 'b', 'EdgeColor', 'b', 'NumBins', 40, 'BinLimits', [0 1], 'Normalization', 'probability', 'FaceAlpha', 0.3);
        histogram(reshape(allData(k).respROIdff0_maxR_sq(2,:,:), 1, []), 'FaceColor', 'r', 'EdgeColor', 'r', 'NumBins', 40, 'BinLimits', [0 1], 'Normalization', 'probability', 'FaceAlpha', 0.3);
        histogram(reshape(allData(k).respROIdff0_maxR_sq(1,:,:), 1, []), 'FaceColor', [255/256,215/256,0], 'EdgeColor', [255/256,215/256,0], 'NumBins', 40, 'BinLimits', [0 1], 'Normalization', 'probability', 'FaceAlpha', 0.3);
        hold off 
        xlabel('max R val')
        ylabel('ROI prop')
        title('hist of maxR')
        suptitle(sprintf('Tad %d (t=%d) maxR by stimtype', allData(k).expnum, k))
        fig_filename = sprintf('Tad %d (t=%d) maxR by stimtype', allData(k).expnum, k)
        saveas(gcf, fig_filename, 'png')
        saveas(gcf, fig_filename, 'epsc2') 
        close;
    end
end

    