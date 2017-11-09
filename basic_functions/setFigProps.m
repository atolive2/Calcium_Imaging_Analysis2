function SetFigProps( gcf )
%SetFigProps sets the figure font sizes and such to be what I want, not
%what Matlab thinks I want.
%   Detailed explanation goes here

set(gcf, 'FontSize', 12);
set(gcf, 'TitleFontSizeMultiplier', 1.3);
colormap(gcf, 'summer')

end