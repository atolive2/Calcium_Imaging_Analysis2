% Examples from GUI exercise guide Simakov 2006

hf=figure;
ha=axes;
x=[-10:0.1:10];
h=line(x,exp(-x.^2));

get(hf,'Type')
get(ha,'Parent')==hf
P=get(ha,'Position');
set(ha,'Position',[P(1) P(2)+P(4)/4 P(3) P(4)/2]);
set(hf,'Color',[0 0 0]);
set(h,'Color',[1 0 0],'Linewidth',5);
hf1=figure; ha1=axes('Parent',hf1);
set(h,'Parent',ha1);

hf=figure;
set(hf,'HandleVisibility','on');
get(hf,'HandleVisibility')
gcf==hf
rootCh=get(0,'Children');
find(rootCh==hf)
set(hf,'HandleVisibility','off');
get(hf,'HandleVisibility')
gcf==hf
rootCh=get(0,'Children');
find(rootCh==hf)

hf=figure;
hL=plot(rand(10,1),'Tag','MyLine');
h=get(hf,'CurrentObject')
get(h,'Tag')
disp(h-hL);

hf=figure('NumberTitle','off');
f=@(h,eD) set(gcbf,'Name',['CallbackObject: ' get(gcbo,'Type')]);
ha=axes('Position',[0.1 0.3 0.8 0.6]);
hL=line([1:10],rand(1,10));
hb=uicontrol('Style','Pushbutton');
hum=uimenu('Label','Test');
set([hf ha hL],'ButtonDownFcn',f);
set([hb hum],'Callback',f);