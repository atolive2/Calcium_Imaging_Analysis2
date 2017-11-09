%% set your directory, either with gui or hand code it:

% first you need to tell it where to get tifs
% if you want gui selection make this bool 1
useGui=1;

% if you want to hand code it, set the above bool to 0 and
% set the path below, do not add the ending slash, I do that later
%imageDirectoryPath = 'F:\Calcium_Imaging_Analysis\_unanalyzed data\20170517 Ca exp 31\concat_trials'
imageDirectoryPath = 'D:\\Torrey_calcium_imaging\20170913 ca exp 42\concat_trials'
%imageDirectoryPath='/Users/cad/Desktop/tAutoTest';
% on non *nix system it will be F:\myImagesEtc

% this actually sets the path
if useGui==1
	imPath=uigetdir();
elseif useGui==0
	imPath=imageDirectoryPath;
end

% lets filter for good measure, looks like your images start with exp
% will work without, but better to be safe.
filterString='exp';
% change the string above if you need.
% below we look at the contents of your directory and get the deets on all
% of the files that meet your sting criterion
filteredFiles = dir([imPath filesep '*' filterString '*']);
% i count the tifs in there
tifCount=numel(filteredFiles);


% make a template
% with this targTemplateStack call, I find the 'middle' tif and then I will
% start looking for template frames int he middle
targTemplateStack=fix(tifCount/2);
tic
for n=targTemplateStack
    mpTifInfo=imfinfo([imPath filesep filteredFiles(n).name]);
    % this counts how many images are in each tif
    imageCount(:,n)=length(mpTifInfo);
    imageBitDepth(:,n)=mpTifInfo(n).BitDepth;
    imageWidths(:,n)=mpTifInfo(n).Width;
    imageHeights(:,n)=mpTifInfo(n).Height;
    fullTifPaths{n}=[imPath filesep filteredFiles(n).name];
    eval(['tempLum=zeros(' num2str(imageCount(:,n)) ',1,''uint' num2str(imageBitDepth(:,n))  ''');']);
    for k=1:imageCount(:,n)
        eval(['templateStack(:,:,' num2str(k) ')=imread(fullTifPaths{' num2str(n) '},''Index'',' num2str(k) ');'])
        eval(['tempLum(' num2str(k) ',1)=mean2(imread(fullTifPaths{' num2str(n) '},''Index'',' num2str(k) '));'])
        if mod(k,200)==0
            disp(['on image ' num2str(k) ' of ' num2str(imageCount(:,n)) ' for stack # ' num2str(n) ' of ' num2str(tifCount) ')'])
        end
    end
end
toc
% the name template stack is temporary


templateFrameNum=200;
tTempStart=fix(numel(tempLum)/2);
halfTempLum=double(tempLum(tTempStart:end,1));
% figure,plot(diff(double(halfTempLum)))
tDiff=diff(halfTempLum);
tMDiff=mean(tDiff);
tMStd=std(tDiff);

ogDiffFrames=numel(tDiff);
candTFrames=find(tDiff<=(tMDiff+tMStd) & tDiff>=(tMDiff-tMStd));
restrictedDiffFrames=numel(candTFrames);

if templateFrameNum>restrictedDiffFrames
    templateFrameNum=restrictedDiffFrames;
end

disp([num2str(ogDiffFrames-restrictedDiffFrames) '/' num2str(ogDiffFrames) ' frames are unfit for template'])
disp(['ill find ' num2str(templateFrameNum) ' frames at random from the rest of the ' num2str(restrictedDiffFrames) ' still frames'])
rng('shuffle')
templateFrames=randperm(templateFrameNum);
templateFrame=im2uint16(mean(templateStack(:,:,templateFrames),3),'Indexed');

figure;
imagesc(templateFrame) 
savefig('templateFrame')
%% now register that stack since you already imported it.
preRegMean=mean(templateStack,3);
regTemp=templateFrame;
subpixelFactor=100;
totalImagesPossible=imageCount(targTemplateStack);
curRegisteredTransformations=zeros(4,totalImagesPossible);
tic
for n=1:totalImagesPossible
    im2reg=templateStack(:,:,n);
    [out1,out2]=dftregistration(fft2(regTemp),fft2(im2reg),subpixelFactor);
    curRegisteredTransformations(:,n)=out1;
    templateStack(:,:,n)=abs(ifft2(out2));
    if mod(n,200)==0
        disp(['registered ' num2str(n) ' images of ' num2str(totalImagesPossible)])
    end
end
toc
postRegMean=mean(templateStack,3);
figure,imshowpair(preRegMean,postRegMean,'diff')

eval(['registeredStack_' num2str(targTemplateStack) '=templateStack;'])
clear templateStack
eval(['registrationTransforms_' num2str(targTemplateStack) '=curRegisteredTransformations;'])
clear curRegisteredTransformations
%%
% now I loop through the other tifs and get th eir name and full path
% note the setdiff thing. I make a list of all the tifs that should be
% there and null the one you already registered.
tic
tifsLeft=tifCount-1;
totalTifVec=1:tifCount
tifsLeft=setdiff(totalTifVec,targTemplateStack)

if numel(tifsLeft)>=1
    for n=tifsLeft
        mpTifInfo=imfinfo([imPath filesep filteredFiles(n).name]);
        % this counts how many images are in each tif
        imageCount(:,n)=length(mpTifInfo);
        imageBitDepth(:,n)=mpTifInfo(n).BitDepth;
        imageWidths(:,n)=mpTifInfo(n).Width;
        imageHeights(:,n)=mpTifInfo(n).Height;
        fullTifPaths{n}=[imPath filesep filteredFiles(n).name];
        eval(['tempStack=zeros(' num2str(imageHeights(:,n)) ',' num2str(imageWidths(:,n)) ',' num2str(imageCount(:,n)) ',''uint' num2str(imageBitDepth(:,n))  ''');']);
        curRegisteredTransformations=zeros(4,imageCount(:,n));
        for k=1:imageCount(:,n)
            im2reg=imread(fullTifPaths{n},'Index',k);
            [out1,out2]=dftregistration(fft2(regTemp),fft2(im2reg),subpixelFactor);
            curRegisteredTransformations(:,k)=out1;
            tempStack(:,:,k)=abs(ifft2(out2));
            if mod(k,200)==0
                disp(['on image ' num2str(k) ' of ' num2str(imageCount(:,n)) ' for stack # ' num2str(n) ' of ' num2str(tifCount) ')'])
            end
        end
        disp(['imported stack # ' num2str(n) ' of ' num2str(tifCount)])
        eval(['registeredStack_' num2str(n) '=tempStack;'])
        clear tempStack
        eval(['registrationTransforms_' num2str(targTemplateStack) '=curRegisteredTransformations;'])
        clear curRegisteredTransformations
    end
    toc
end
 

%% concatenate stacks I haven't fully tested this, so you may need to tweak
if tifCount>1
catStack=registeredStack_1;
clear registeredStack_1
    for n=2:tifCount
        eval(['catStack=cat(3,catStack,registeredStack_' num2str(n) ');']);
        eval(['clear registeredStack_' num2str(n)])
    end
    %eval(['clear registeredStack_' num2str(n)])
end

%% If more than 6 blocks, split catStack into 2 for roiInspector
catStack1 = catStack(:,:,1:10000);
catStack2 = catStack(:,:,10001:end);
%then delete catStack
%then run roiInspector extraction on catStack1 and rename somaticF to
%somaticF1
%then run roiInspector extraction on catStack2 and rename somaticF to
%somaticF2

%% exp 42 problem with registration testing

catStack_sm = cat(3, catStack(:,:,1:3840), catStack(:,:,5761:19200));

% concatenate specific registered stacks
keep = [1 4 5 6 7];
catStack=registeredStack_1;
%clear registeredStack_1
    for n=2:length(keep)
        eval(['catStack=cat(3,catStack,registeredStack_' num2str(keep(n)) ');'])
        %sprintf('registeredStack_%s', num2str(keep(n)))
        %eval(['clear registeredStack_' num2str(n)])
    end

% average the average images
for k = 1:length(keep)
    %'avg_keep(:,:,k) = meanProj_registeredStack_' num2str(keep(n)) ');']);

    eval(['avg_keep(:,:,k) = meanProj_registeredStack_' num2str(keep(n)) ');'])
    
end

