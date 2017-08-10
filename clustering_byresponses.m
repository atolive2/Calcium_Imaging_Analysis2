%% Probability of response as a way to define networks

% Asking these questions to see if there are cells that fire together (fire
% together wire together adage)
    %What's the probability of responding given that this cell responds? 
    %Are there any cells with the same.activation pattern(eg respond in all the same trials) 
    
%% Define "Responded" by peak val > 0.015

for t = 1:length(tadcluster)
    for i = 1:size(tadcluster{1,t}.peak_bytrial,1)
        for j = 1:size(tadcluster{1,t}.peak_bytrial,2)
            if tadcluster{1,t}.peak_bytrial(i,j) > 0.015
                tadcluster{1,t}.response_bypeak(i,j) = 1
            else
                tadcluster{1,t}.response_bypeak(i,j) = 0
            end
        end
    end
end

% Are there any identical rows?

for t = 1:length(tadcluster)
    tadcluster{1,t}.unique_rows = unique(tadcluster{1,t}.response_bypeak, 'rows')
end

for t = 1:length(tadcluster)
    if size(tadcluster{1,t}.response_bypeak) == size(tadcluster{1,t}.unique_rows)
        any_unique(t) = 1
    end
end

%no, there are not any unique rows in any experiments (2,3,4 have rows of
%zeros)

% Are there any ROI pairs with at least 50% agreement? 

for t = 1:length(tadcluster)
for i = 1:size(tadcluster{1,t}.response_bypeak,1)
    for j = 1:size(tadcluster{1,t}.response_bypeak,1)
        for k = 1:size(tadcluster{1,t}.response_bypeak,2)
            if tadcluster{1,t}.response_bypeak(i,k) == tadcluster{1,t}.response_bypeak(j,k)
                tadcluster{1,t}.same_response(i,k) = 1;
            else
                tadcluster{1,t}.same_response(i,k) = 0;
            end
        end
        tadcluster{1,t}.sum_same_response(i, j) = sum(tadcluster{1,t}.same_response(i,:));
    end
end
end

for t = 1:length(tadcluster)
    tadcluster{1,t}.prop_same_response = tadcluster{1,t}.sum_same_response ./ size(tadcluster{1,t}.response_bypeak,2);
end

% plot prop same response as a heat map
for t = 1:length(tadcluster)
    figure;
    colormap('hot')
    colormap(flipud(colormap))
    imagesc(tadcluster{1,t}.prop_same_response)
    colorbar
    title(sprintf('tad %d prop agreement between ROIs', t))
    fig_filename = sprintf('tad %d prop agreement between ROIs', t)
    saveas(gcf,fig_filename,'png');
    close;
end

% Proportion in ranges 0-.25, .25-.5, .5-.75, .75-0.999
for t = 1:length(tadcluster)
    prop_rng(t,1) = length(find(tadcluster{1,t}.prop_same_response < 0.25));
    prop_rng(t,2) = length(find( 0.25 <= tadcluster{1,t}.prop_same_response < 0.5));
    prop_rng(t,3) = length(find( 0.5 <= tadcluster{1,t}.prop_same_response < 0.75));
    prop_rng(t,3) = length(find( 0.75 <= tadcluster{1,t}.prop_same_response <= 1));
    prop_rng(t,4) = size(tadcluster{1,t}.prop_same_response,1);
    prop_rng(t,5) = prop_rng(t,4)^2 - prop_rng(t,4);
end

for t = 1:size(prop_rng,1)
    prop_toplot(t,:) = prop_rng(t,1:4)./(prop_rng(t,4)^2)
end

bar(prop_toplot, 'stacked')
legend({'0-.25', '.25-.5', '.5-.75', '.75-1'}, 'Orientation', 'horizontal')