%% What is special about the the highly correlated cells
% Use cluster_data_respROI with 0=not highly corr, 1=highly corr

% stat test all vars separated by highly/not highly corr

highly_corr = find(cluster_data_respROI(:,22))
not_highly_corr = find(~cluster_data_respROI(:,22))

for v = 3:(size(cluster_data_respROI,2)-1)
     [P_KW(v), TBL_KW{v}, stats_KW{v}] = kruskalwallis(cluster_data_respROI(:,v), cluster_data_respROI(:,22));
end

sig_vars = find(P_KW < 0.05)
labels = {'peak MS', 'peak V', 'peak M', 'peak NS', 'peak SD MS', 'peak SD V', 'peak SD M', 'peak SD NS', 'onset time MS', 'onset time V', 'onset time M', 'onset time NS', ...
    'onset time SD MS', 'onset time SD V', 'onset time SD M', 'onset time SD NS', ...
    'MSEnh peak', 'Uni bias num resp'}

for v = 1:length(sig_vars)
    figure;
    subplot(1,2,1)
    hist(cluster_data_respROI(highly_corr, sig_vars(v)))
    title('highly correlated')
    
    subplot(1,2,2)
    hist(cluster_data_respROI(not_highly_corr, sig_vars(v)))
    title('not highly correlated')
    
    suptitle(sprintf('Var %s sig P_{KW} between highly/not highly corr', labels{sig_vars(v)-2}))
    fig_filename = (sprintf('Var %s sig P_{KW} between highly and not highly corr', labels{sig_vars(v)-2}))
    saveas(gcf, fig_filename, 'png')
    close;
end
