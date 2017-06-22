%% Data Viewer GUI for Carlos

%% Section 1: tab through all experiments, see what they are, view individual ROIs

h = figure; 
get(h) % get figure properties

% Create the navigator from exp to exp

% create an editable textbox object for typing in a particular exp
edit_box_h = uicontrol('style','edit',...
                    'units', 'normalized',...
                    'position', [0.15 0.95 0.05 0.05]);
% create a button to go to previous exp
but_h1 = uicontrol('style', 'pushbutton',...
                    'string', 'next exp',...
                    'units', 'normalized',...
                    'position', [0.2 0.95 0.1 0.05])
                    %'callback', {@eg_fun,edit_box_h, ellipse_h });
% create a button to go to next exp
but_h2 = uicontrol('style', 'pushbutton',...
                    'string', 'prev exp',...
                    'units', 'normalized',...
                    'position', [0.05 0.95 0.1 0.05])
                    %'callback', {@eg_fun,edit_box_h, ellipse_h });                    
                    
% Create and Position Graphs 
%[in from left side, bottom x axis loc, end of graph right side, end of graph top]

% top left/tectum scatter of primary modality
ha = axes;
set(ha, 'Position', [0.05 0.55 0.3 0.3])
% top center/tectum scatter of time to peak gradient
hb = axes;
set(hb, 'Position', [0.4 0.55 0.3 0.3])
% top right/tectum scatter of variability of time to peak
hc=axes;
set(hc,'Position',[0.75 0.55 0.3 0.3])

% bottom left/stacked bar of primary modality split
hd=axes;
set(hd,'Position',[0.05 0.05 0.1 0.35])
% bottom center/ all traces for a given cell (by stim modality)
he = axes;
set(he, 'Position', [0.2 0.05 0.4 0.35])
% bottom right/1 point of each trace by stim modality
hf = axes;
set(hf, 'Position', [0.68 0.05 0.3 0.35])

% Create Graph Labels (as text box annotations)
an_a = annotation('textbox',[0.05 0.82 0.3 0.1],'String','Primary Modality','FitBoxToText','on');
an_b = annotation('textbox',[0.4 0.82 0.3 0.1],'String','Time to Peak','FitBoxToText','on');
an_c = annotation('textbox',[0.75 0.82 0.3 0.1],'String','var of time to peak','FitBoxToText','on');
an_d = annotation('textbox',[0.05 0.38 0.3 0.1],'String','Prop PM','FitBoxToText','on');
an_e = annotation('textbox',[0.23 0.38 0.3 0.1],'String','All Sweeps','FitBoxToText','on');
an_f = annotation('textbox',[0.7 0.38 0.3 0.1],'String','Peak by Sweep','FitBoxToText','on');

%% Graph contents
set(ha, 'UserData', tadpole)

%% Collect data into easy to use, obvious vars in each tadpole


% get roi centers in a matrix
for t = 1:length(tadpole)
    for r = 1:length(tadpole{1,t}.somaticROICenters)
        tadpole{1,t}.ROIlocs_dbl(r,:) = tadpole{1,t}.somaticROICenters{1,r}(1).Centroid;
    end
end
for t = 1:length(tadpole)
    tadpole{1,t}.ROIlocs_dbl(:,2) = 250 + (250 - tadpole{1,t}.ROIlocs_dbl(:,2));
end

% Get colors for ROIs - primary modality
for i = 1:

% get colors for ROIs - time to peak (use Jet)

% get colors for variability of time to peak (use Jet)





%%%%%%%%%%%%%
ha.Children = scatter(tadpole{1,t}.ROIlocs_dbl(:,1), tadpole{1,t}.ROIlocs_dbl(:,2), 'k')

hb.UserData = scatter(tadpole{1,t}.ROIlocs_dbl(:,1), tadpole{1,t}.ROIlocs_dbl(:,2), 'r')