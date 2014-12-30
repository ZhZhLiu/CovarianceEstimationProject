function [] = ExportRhoEvolutiton(rhoChain,fileName)

nSim = length(rhoChain);

handle = figure(4);

plot(rhoChain);
str = sprintf('\\rho^2');
title(str);
ylim([0, 1]);
xlim([0,nSim]);

nameFile_EPS = sprintf('RhoTestFigures/Eps/%s.eps', fileName);
nameFile_FIG = sprintf('RhoTestFigures/Fig/%s.fig', fileName);
nameFile_JPG = sprintf('RhoTestFigures/Jpg/%s.jpg', fileName);

set(gcf,'PaperPositionMode','auto')
set(handle, 'Position', [100, 100, 1400, 900]);
print(handle, '-depsc', nameFile_EPS);
print(handle, '-djpeg', nameFile_JPG);
savefig(handle,nameFile_FIG);


end

