%% Automation of importer registration
% This code uses parts of Chris Deister's GUI to automate the import,
% registration and saving of each .tiff stack individually (rather than
% manually creating a large stack and running it all at once which creates 
% RAM problems above 5 trial blocks). 

%% Import one multi-page tif that you think will give you the best template, and make that template. 
% Use the GUI for this part. 

% Create a matrix with the get path information in it.
Path = ['F:/Calcium_Imaging_Analysis/_unanalyzed data/20170517 Ca exp 31/concat_trials/']


%% Import each multi-page tif one by one, and keep the same template.
% start looping through each stack of tiffs.
exp = 31
blocks = [2 3 4 5 6 7 8]
b=1
importPath = sprintf([Path, 'exp%dblock%d'], exp, blocks(b)) 

% import image stack (I kept everything from the importer, but I only need
% multipage tiffs)

% from lines 67-204 of importer.m
% --- Executes on button press in importButton.
function importButton_Callback(hObject, eventdata, handles)
% hObject    handle to importButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

mPF=get(handles.multiPageFlag, 'Value');
pImport=get(handles.parallelizeImportToggle,'Value');

% Check to see if the user imported something already and/or wants to
% import a multi-page Tif. 
g=evalin('base','exist(''importPath'')');
disp('importing images ...')

% User has set a path, but doesn't want multi-page tif.
if g==1 && mPF==0
    imPath=evalin('base','importPath');
    firstIm=str2num(get(handles.firstImageEntry,'string'));
    endIm=str2num(get(handles.endImageEntry,'string'));
    
% User has not set a path, and doesn't want multi-page tif.    
elseif g==0 && mPF==0
    imPath=uigetdir;
    firstIm=str2num(get(handles.firstImageEntry,'string'));
    endIm=str2num(get(handles.endImageEntry,'string'));
    
% User has set a path, but wants multi-page tif.    
elseif g==1 && mPF==1
    imPath=evalin('base','importPath');
    tifFile=evalin('base','tifFile');
    mpTifInfo=evalin('base','mpTifInfo');
    firstIm=str2num(get(handles.firstImageEntry,'string'));
    endIm=str2num(get(handles.endImageEntry,'string'));
    
% User has set not path, and does want multi-page tif.    
elseif g==0 && mPF==1
    [tifFile,imPath]=uigetfile('*.*','Select your tif file');
    mpTifInfo=imfinfo([imPath tifFile]);
    imageCount=length(mpTifInfo);
    firstIm=str2num(get(handles.firstImageEntry,'string'));
    endIm=str2num(get(handles.endImageEntry,'string'));
    assignin('base','mpTifInfo',mpTifInfo);
    assignin('base','tifFile',tifFile);
    assignin('base','imPath',imPath);
end

% This loads a file list that has characters that match the filter string.
% It should detect the bit depth and dimensions.
if mPF==0
    filterString={get(handles.fileFilterString,'String')};
    filteredFiles = dir([imPath filesep '*' filterString{1} '*']);
    filteredFiles=resortImageFileMap(filteredFiles);
    assignin('base','filteredFiles',filteredFiles)
    importCount=(endIm-firstIm)+1;
    disp(imPath)
    disp(filteredFiles(1,1).name)
    canaryImport=imread([imPath filesep filteredFiles(1,1).name]);
    imageSize=size(canaryImport);
    canaryInfo=whos('canaryImport');
    bitD=canaryInfo.class;
    assignin('base','bitDebug',bitD); % debug 
    importedImages=zeros(imageSize(1),imageSize(2),importCount,bitD);
    if strcmp(bitD,'uint16')==1
        imType='uint16';
    elseif strcmp(bitD,'uint32')==1
        imType='uint32';
    elseif strcmp(bitD,'uint8')==1
        imType='uint8';
    else
        imType='Double';
    end
    disp(imType)
 
    
    tic
    if pImport==1
        tempFiltFiles=filteredFiles(firstIm:endIm,1);
        parfor n=1:importCount;
            importedImages(:,:,n)=imread([imPath filesep tempFiltFiles(n,1).name]);
        end
    elseif pImport==0
        for n=firstIm:endIm;
            importedImages(:,:,(n-firstIm)+1)=imread([imPath filesep filteredFiles(n,1).name]);
        end
    end
    iT=toc;
    
    if bitD==16
        assignin('base',['importedStack_' filterString{1}],uint16(importedImages));
    elseif bitD==8
        assignin('base',['importedStack_' filterString{1}],uint8(importedImages));
    elseif bitD==32
        assignin('base',['importedStack_' filterString{1}],uint32(importedImages));
    else
        assignin('base',['importedStack_' filterString{1}],double(importedImages));
    end
    vars = evalin('base','who');
    set(handles.workspaceVarBox,'String',vars)
    
else  % The user wants multi-page tif. This import is a bit different.
    bitD=mpTifInfo(1).BitDepth;
    mImage=mpTifInfo(1).Width;
    nImage=mpTifInfo(1).Height;
    NumberImages=length(mpTifInfo);
    if bitD==16
        imType='uint16';
    elseif bitD==32
        imType='uint32';
        % why are you using 32 bit images? I'm curious shoot me an email please.
    elseif bitD==8
        imType='uint8';
    else
        imType='Double';
    end
 
    importedStack=zeros(nImage,mImage,NumberImages,imType);
    tic
    if pImport==1
        parfor i=1:NumberImages
            importedStack(:,:,i)=imread([imPath tifFile],'Index',i);
        end
    elseif pImport==0
        for i=1:NumberImages
            importedStack(:,:,i)=imread([imPath tifFile],'Index',i);
        end
    end
    assignin('base','importedStack',importedStack)
    assignin('base','importedBitDepth',bitD)
    iT=toc;
    
    % update var box
    vars = evalin('base','who');
    set(handles.workspaceVarBox,'String',vars)
end

disp(['*** done with import, which took ' num2str(iT) ' seconds'])
% Update handles structure
guidata(hObject, handles);


% Set the newly imported stack as the registration stack.

% from lines 454-475 of importer.m
% --- Executes on button press in setRegStackButton.
function setRegStackButton_Callback(hObject, eventdata, handles)
% hObject    handle to setRegStackButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
selections = get(handles.workspaceVarBox,'String');
selectionsIndex = get(handles.workspaceVarBox,'Value');
if numel(selectionsIndex)>1
    disp('you can only set one primary stack to register')
elseif numel(selectionsIndex)==0
    disp('you must select a stack to register')
else
    assignin('base','stackToRegister',selections{selectionsIndex});
end

vars = evalin('base','who');
set(handles.workspaceVarBox,'String',vars)
    


% Update handles structure
guidata(hObject, handles);


% Change registration stack. You may need to append each stack with a number.

% Register stack.
% from lines 481-526 of importer.m
% --- Executes on button press in registerButton.
function registerButton_Callback(hObject, eventdata, handles)
% hObject    handle to registerButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% get the string for the stack you want to register; this will be a string
regStackString=evalin('base','stackToRegister');

% todo: allow user to crop what they want
% the way I do the registration rotates the image, here I offset that


regTemp=evalin('base','regTemplate');
rStack=evalin('base',regStackString);
subpixelFactor=100;
totalImagesPossible=size(rStack,3);

% pre-allocate, because ... matlab ...

registeredImages=zeros(size(rStack,1),size(rStack,2),totalImagesPossible,'uint16');
registeredTransformations=zeros(4,totalImagesPossible);

disp('registration started ...')

tic
regTempC=regTemp;
    parfor n=1:totalImagesPossible,
        imReg=rStack(:,:,n);
        [out1,out2]=dftregistration(fft2(regTempC),fft2(imReg),subpixelFactor);
        registeredTransformations(:,n)=out1;
        registeredImages(:,:,n)=abs(ifft2(out2));
        %registeredImages(:,:,n)=imrotate(abs(ifft2(out2)),180);
    end
t=toc;
disp(['done with registration. it took ' num2str(t) ' seconds'])
assignin('base',[regStackString '_registered'],uint16(registeredImages))
assignin('base','registeredTransforms',registeredTransformations)


% update var box
vars = evalin('base','who');
set(handles.workspaceVarBox,'String',vars)

% Update handles structure
guidata(hObject, handles);

% Delete the unregistered imported stack.

% Import the next multi-page tif, set it as reg stack and register.

% Delete the unregistered

% Concatenate registered stacks 