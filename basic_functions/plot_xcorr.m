function plot_xcorr( data, color )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
figure;
colormap(color)
colormap(flipud(colormap))


imagesc(abs(data))
colorbar
ylabel('ROI', 'fontsize', 30)
xlabel('ROI', 'fontsize', 30)
set(gca, 'fontsize', 20)
end

