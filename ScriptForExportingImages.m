%insert this script figname after whichever figure you want to export from
%MATLAB to tiff 
% Change the figname accordingly, then use Image Magick for the last step
% from the command prompt

% Before calling this script for a figure, type the following for each
% figure:
% figname = 'figure2'; % change the figure name accordingly

% Use the pdfs in the manuscript and the tiff files for the journal
fig= gcf;

% Set the desired output size in inches (e.g., 6 inches wide full text width)
% widthInches = 20; heightInches = 10; for 2x4 grid 
% widthInches = 15; heightInches = 10; for 2x3 grid or 1x3
% as height scales by export graphics
% widthInches = 10; heightInches = 10; for 2x2 grid 
%widthInches = 10;
%heightInches = 10; % Adjust based on your figure's aspect ratio


% Get current size in inches
set(fig, 'Units', 'Inches');
figPos = get(fig, 'Position');


% Set figure size
set(fig, 'Position', [1, 1, widthInches, heightInches]);
set(fig, 'PaperUnits', 'Inches');
set(fig, 'PaperSize', [widthInches, heightInches]);
figPos = get(fig, 'Position');
set(fig, 'PaperPositionMode', 'manual');
fig = gcf;


% Apply to all subplots
ax = findall(gcf, 'Type', 'axes');
for k = 1:length(ax)
    box(ax(k), 'on'); % or 'off' depending on your preference
end

exportgraphics(fig,[figname, '.tiff'],'Resolution',600)
exportgraphics(fig,[figname, '.pdf'],'Resolution',600)

% do the following in the command prompt with image magick installed
%magick figureS4big.tiff -resize 34.48% -density 600 figureS4.tiff

% Desired output = 3 inches (approximately 1 column)
% 2panel figure 8.7x6.64, scale 3/8.7 = 34.48%