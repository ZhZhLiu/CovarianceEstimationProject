function [] = ExportRhoEvolutiton(rhoChain,fileName, rootPath)

nSim = length(rhoChain);

% Handling the directory and folder

epsPath = strcat(rootPath, 'MatrixSurfacePlot/EPS/');
figPath = strcat(rootPath, 'MatrixSurfacePlot/FIG/');
jpgPath = strcat(rootPath, 'MatrixSurfacePlot/JPG/');

nameFile_EPS = strcat(epsPath, sprintf('%s.eps', fileName));
nameFile_FIG = strcat(figPath, sprintf('%s.fig', fileName));
nameFile_JPG = strcat(jpgPath, sprintf('%s.jpg', fileName));

if exist(epsPath,'dir') == 0
    mkdir(epsPath);
end

if exist(figPath,'dir') == 0
    mkdir(figPath);
end

if exist(jpgPath,'dir') == 0
    mkdir(jpgPath);
end

% Plot the graph

handle = figure(4);

plot(rhoChain);
str = sprintf('\\rho^2');
title(str);
ylim([0, 1]);
xlim([0,nSim]);

% Export the graph

set(gcf,'PaperPositionMode','auto')
set(handle, 'Position', [100, 100, 1400, 900]);
print(handle, '-depsc', nameFile_EPS);
print(handle, '-djpeg', nameFile_JPG);
savefig(handle,nameFile_FIG);

end

