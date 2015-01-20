function [] = ExportSigma2Evolutiton(sig2Chain,fileName, rootPath)

[p,nSim] = size(sig2Chain);

handle = figure(1);

% Handling the directory and folder

epsPath = strcat(rootPath, 'VarianceSimPaths/EPS/');
figPath = strcat(rootPath, 'VarianceSimPaths/FIG/');
jpgPath = strcat(rootPath, 'VarianceSimPaths/JPG/');

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

% Plotting the graph

plotRowIndices = floor(linspace(2, p, 9));
maxRange =  max(max(sig2Chain(plotRowIndices, :)));
minRange =  min(min(sig2Chain(plotRowIndices, :)));

for i = 1:length(plotRowIndices)
    subplot(3,3,i);
    rowIndex = plotRowIndices(i);
    plot(sig2Chain(rowIndex, :));
    str = sprintf('\\sigma^2 with k = %d',rowIndex);
    title(str);
    ylim([minRange, maxRange]);
    xlim([0,nSim]);
end

% Save graph to file
set(gcf,'PaperPositionMode','auto')
set(handle, 'Position', [100, 100, 1400, 900]);
print(handle, '-depsc', nameFile_EPS);
print(handle, '-djpeg', nameFile_JPG);
savefig(handle,nameFile_FIG);

end

