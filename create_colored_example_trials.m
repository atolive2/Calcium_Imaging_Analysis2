% Use grs2rbg function to convert a .tiff stack to RGB pseudocolor
% For Carlos to send to article writer at Nature ?

% http://www.mathworks.com/matlabcentral/fileexchange/13312-grayscale-to-rgb-converter

% the function grs2rbg takes a single image - not a tif stack. Therefore,
% loop through an imported tiff file. The code below assumes you have used
% the ImageAnalysis_GUI from Chris Deister to import a .tiff stack into
% importedStack. THis code takes a long time to run, so stick with 1 trial
% at a time. 

% built in colormaps: https://www.mathworks.com/matlabcentral/mlc-downloads/downloads/submissions/35242/versions/3/screenshot.jpg

% original version - with parula heatmap.
cmap = colormap(parula);
for i = 1:size(importedStack, 3)
    Res(:,:,:,i) = grs2rgb(importedStack(:,:,i), cmap); 
end

image(Res) % to inspect the coloration of an image

% try different colormaps
cmap = colormap(jet);
for i = 1:size(importedStack, 3)
    Res_lines(:,:,:,i) = grs2rgb(importedStack(:,:,i), cmap); 
end

image(Res_lines(:,:,:,1))

% write RGB matrix with dims img_size x img_size x 3 x num_images to file
% (creates a new .tiff). 
imwrite(Res_lines(:,:,:,1), 'e5b6_st1281_jet.tiff')
for i = 2:size(Res_lines, 4)
    imwrite(Res_lines(:,:,:,i),'e5b6_st1281_jet.tiff','WriteMode','append');
end

% then open this newly created .tiff in Fiji (is just ImageJ) to convert to
% .avi with desired framerate.
% it seems like the matlab function videoWriter might do this part, but I
% can't get it to work. 