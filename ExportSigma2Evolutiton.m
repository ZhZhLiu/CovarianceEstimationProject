function [] = ExportSigma2Evolutiton(sig2Chain,fileName)

[p,nSim] = size(sig2Chain);

handle = figure(1);

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

nameFile_EPS = sprintf('VarianceTestFigures/Eps/%s.eps', fileName);
nameFile_FIG = sprintf('VarianceTestFigures/Fig/%s.fig', fileName);
nameFile_JPG = sprintf('VarianceTestFigures/Jpg/%s.jpg', fileName);

set(gcf,'PaperPositionMode','auto')
set(handle, 'Position', [100, 100, 1400, 900]);
print(handle, '-depsc', nameFile_EPS);
print(handle, '-djpeg', nameFile_JPG);
savefig(handle,nameFile_FIG);

end

