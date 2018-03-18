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

    
%% Next step (For AAAS 2018 poster)
% create paired images for sample ROIs
% convientiently, all ROIs of tad 34 respond, so respROI = actual ROI

t = 17 %exp 34
roi = 2 % a highcorr ROI

for r = 1:length(tadcluster{1,t}.resp_ROIs)
    if r == roi
        continue
    else
        figure;
    hold on
    
    plot(tadcluster{1,t}.dff0_multi(roi,:), 'Color', 'k', 'LineWidth', 2)
    plot(tadcluster{1,t}.dff0_multi(r,:)+0.3, 'Color', H(roi,r,:), 'LineWidth', 2)

    hold off
    px = round(pdist2(tadcluster{1,t}.ROIcenters(roi,:), tadcluster{1,t}.ROIcenters(r,:), 'euclidean'), 0)
    xlim([636 2385])
    ylim([-1 1])
    ax = gca;
    ax.Visible = 'off'
%     title(sprintf('Tad %d ROI %d vs ROI %d (%d px apart)', t, roi, r, px), 'fontsize', 30)
%     ylabel('\DeltaF/F_{0}', 'fontsize', 30)
%     xlabel('time (frames)', 'fontsize', 30)
%     set(gca, 'fontsize', 30)
    fig = gcf;
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 8 3];
    fig_filename = sprintf('Tad %d ROI %d vs ROI %d (%d px apart) no axes', tadcluster{1,t}.expnum, roi, r, px)
    saveas(gcf, fig_filename, 'png')
    
    close;
        
    end
end

%% Get specific trials for specific ROIs

%% For methods examples
% ROI 3
% plot all MS trials on top of each other
MSstims = [1 5 8 9]
Vstims = [2 6]
Mstims = [3 7]

roi = 3
figure;
hold on
for s = 1:length(tadcluster{1,t}.stimorder)
    if ismember(tadcluster{1,t}.stimorder(s), MSstims)
        if max(tadcluster{1,t}.df_f0{roi,s}) < 1
            if min(tadcluster{1,t}.df_f0{roi,s}) > -0.1
                plot(tadcluster{1,t}.df_f0{roi,s}, 'LineWidth', 1)
            end
        end
    end
end
hold off
    ax=gca;
    xsize = length(tadcluster{1,t}.df_f0{roi,s});
    ax.XTick = [0, xsize/7, (xsize/7)*2, (xsize/7)*3, (xsize/7)*4, (xsize/7)*5, (xsize/7)*6, (xsize/7)*7];
    ax.XTickLabel = {'0','1', '2', '3', '4', '5', '6', '7'};
    ylabel('\DeltaF/F_{0}', 'fontsize', 30)
    xlabel('time (sec)', 'fontsize', 30)
    set(gca, 'fontsize', 20)
    fig_filename = sprintf('Tad 34 ROI %d all good MS trials df_f0)', roi)
    saveas(gcf, fig_filename, 'png')
    
% ROI 63
roi = 63
figure;
hold on
for s = 1:length(tadcluster{1,t}.stimorder)
    if ismember(tadcluster{1,t}.stimorder(s), MSstims)
        if max(tadcluster{1,t}.df_f0{roi,s}) < 1
            if min(tadcluster{1,t}.df_f0{roi,s}) > -0.1
                plot(tadcluster{1,t}.df_f0{roi,s}, 'LineWidth', 1)
            end
        end
    end
end
hold off
    ax=gca;
    xsize = length(tadcluster{1,t}.df_f0{roi,s});
    ax.XTick = [0, xsize/7, (xsize/7)*2, (xsize/7)*3, (xsize/7)*4, (xsize/7)*5, (xsize/7)*6, (xsize/7)*7];
    ax.XTickLabel = {'0','1', '2', '3', '4', '5', '6', '7'};
    ylabel('\DeltaF/F_{0}', 'fontsize', 30)
    xlabel('time (sec)', 'fontsize', 30)
    set(gca, 'fontsize', 20)
    fig_filename = sprintf('Tad 34 ROI %d all good MS trials df_f0)', roi)
    saveas(gcf, fig_filename, 'png')
    
%% ROIs for correlation examples

%% First, get sorted maxR matrix
% sort by number of highly correlated cells
[B, I] = sort(tadcluster{1,t}.highcorr_numROIs_MS, 'descend')
for dim1 = 1:length(I)
    for dim2 = 1:length(I)
        xcorr_sorted(dim1, dim2) = tadcluster{1,t}.respROIdff0_maxR_sq_MS(I(dim1), I(dim2));
    end
end

figure;
colormap('hot')
colormap(flipud(colormap))
imagesc(xcorr_sorted)
colorbar

ylabel('ROI', 'fontsize', 30)
xlabel('ROI', 'fontsize', 30)
set(gca, 'fontsize', 20)
fig_filename = sprintf('Tad 34 maxR sorted')
saveas(gcf, fig_filename, 'png')
    
    
%% Quantify how many high corr ROIs (for AAAS poster)
% get proportion of cells with at least 25% high corrs per tadpole in each
% modality and scatterplot

% USING DIFFERENT DATA FILE ('allData_workspace_20171109.mat')
for t = 1:length(allData)
    if isfield(allData{1,t}, 'highcorr_numROIs_MS')
        highcorr_quent(1,t) = length(find(allData{1,t}.highcorr_numROIs_MS > (length(allData{1,t}.highcorr_numROIs_MS) / 4)));
        highcorr_quent(2,t) = length(find(allData{1,t}.highcorr_numROIs_V > (length(allData{1,t}.highcorr_numROIs_V) / 4)));
        highcorr_quent(3,t) = length(find(allData{1,t}.highcorr_numROIs_M > (length(allData{1,t}.highcorr_numROIs_M) / 4)));
        highcorr_quant(4,t) = allData{1,t}.stage;
        highcorr_quant(5,t) = length(allData{1,t}.highcorr_numROIs_M);
        
    end
end

% get proportion from raw
highcorr_quantProp(1,:) = highcorr_quant(1,:) ./ highcorr_quant(5,:)

%% Exp 34 ROIs quantified

% make a bar graph of the count of highcorr ROIs for each type
data = [49 33 52 23 34 26 21]
bar(data)
ylabel('ROI count', 'fontsize', 30)
xlabel('Modality', 'fontsize', 30)
set(gca, 'fontsize', 20)
fig_filename = 'Tad 34 bar graph of highcorr ROIs by mod'
saveas(gcf, fig_filename, 'png')