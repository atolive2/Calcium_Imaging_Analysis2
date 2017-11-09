%% Show examples for xcorr
% using exp 34 (t = 17)

t=17

%% Plot the peak sizes for each respROI

%get peak data for each ROI
x_vals = [];
y_vals = [];
for r = 1:length(tadcluster{1,t}.resp_ROIs)
    num_trials = length(tadcluster{1,t}.stimorder);
    x_vals = [x_vals, 2*r*ones(1,num_trials)];
    y_vals = [y_vals, tadcluster{1,t}.peak_bytrial(r,:)];
end

%assign colors based on trial type
    %1 = multisensory high M / high V
    %2 = visual high/crash
    %3 = mechanosensory high
    %4 = no stimulus
    
    %5 = multisensory high M / low V
    %6 = visual low/scrambled crash
    %7 = mechanosensory low
    
    %8 = multisensory low M / high V
    %9 = multisensory low M / low V 

tts = unique(tadcluster{1,t}.stimorder)
MSstims = [1 5 8 9]
Vstims = [2 6]
Mstims = [3 7]
for s = 1:length(tadcluster{1,t}.stimorder)
    if ismember(tadcluster{1,t}.stimorder(s), MSstims) 
        colors(:,s) = [1 0 1] %multi = magenta
    elseif ismember(tadcluster{1,t}.stimorder(s), Vstims) 
        colors(:,s) = [0 0 1] %vis = blue
    elseif ismember(tadcluster{1,t}.stimorder(s), Mstims) 
        colors(:,s) = [1 0 0] %mech = red
    else
        colors(:,s) = [0.5 0.5 0.5] %no stim = gray
    end
end

colors_all = repmat(colors, 1, length(tadcluster{1,t}.resp_ROIs));

%make scatterplot
scatter(x_vals, y_vals, [], colors_all')
xlabel('ROI')
ylabel('peak \DeltaF/F_{0}')
ax = gca;
ax.XTick([1:20:(2*length(tadcluster{1,t}.resp_ROIs))])
ax.XTickLabel(mat2cell(1:10:length(tadcluster{1,t}.resp_ROIs)))

%% Plot total number of responses by ROI
bar(tadcluster{1,t}.sum_repsonses)
ylabel('number of responses')
xlabel('ROI')
axis([-inf inf 0 length(tadcluster{1,t}.stimorder)])

%% Plot example traces

% specifically, plot the whole multisensory timecourse for each ROI, 1 plot
% per ROI has: 1) bold black = ROI of interest 2) each other ROI's trace
% colored by Rval
% offset traces in Y by 0.2

% Define colors based on Rvals (H contains this info)
% Say this is the given matrix:
G = tadcluster{1,t}.respROIdff0_maxR_sq_MS;
% Use IMAGESC to plot G.
colormap(hot) % You realize this affects final image (H)?
colormap(flipud(colormap))
imagesc(G);
title('IMAGESC (MxN)')
% Now make an RGB image that matches display from IMAGESC:
C = colormap;  % Get the figure's colormap.
L = size(C,1);
% Scale the matrix to the range of the map.
Gs = round(interp1(linspace(min(G(:)),max(G(:)),L),1:L,G));
H = reshape(C(Gs,:),[size(Gs) 3]); % Make RGB image from scaled.
figure;
image(H)  % Does this image match the other one?
title('IMAGE (MxNx3)')

% loop through each ROI, plotting
eliminate = [1 3 7 12 13 29 33 35 41 50 60 73];
for r = 1:length(tadcluster{1,t}.resp_ROIs)
    figure;
    hold on
    for i = 1:length(tadcluster{1,t}.resp_ROIs)
        if i == r
            plot(tadcluster{1,t}.dff0_multi(i,:)+(i*0.3), 'Color', H(r,i,:), 'LineWidth', 2)
        elseif ismember(i, eliminate)
            continue
        else 
            plot(tadcluster{1,t}.dff0_multi(i,:)+(i*0.3), 'Color', H(r,i,:))
            %fprintf(num2str(i))
            %pause
        end
    end
    hold off
    title(sprintf('Tad %d xcorr maxR with ROI %d', t, r))
    ylabel('\DeltaF/F_{0}')
    xlabel('time (frames)')
    fig = gcf;
    fig.PaperUnits = 'inches'
    fig.PaperPosition = [0 0 8 3]
    fig_filename = sprintf('Tad %d xcorr maxR with ROI %d', tadcluster{1,t}.expnum, r)
    saveas(gcf, fig_filename, 'png')
    close;
end

    
    

