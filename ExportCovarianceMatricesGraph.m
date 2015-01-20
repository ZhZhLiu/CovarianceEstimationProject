function [] = ExportCovarianceMatricesGraph(SIG_TRUE, COV_EST, SIGINV, COVINV_EST, data, lambda, fileName, rootPath)

handle = figure(2);


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

% Plotting Graph

subplot(2,3,1)
imagesc(SIG_TRUE);
title('True Covariance')
axis square

subplot(2,3,2)
imagesc(cov(data));
title('Sample Covariance')
axis square

subplot(2,3,3)
imagesc(COV_EST);
title(sprintf('Regularized Covariance \\lambda=%4.2f', lambda))
axis square

subplot(2,3,4)
imagesc(SIGINV);
title('True Covariance Inverse')
axis square

subplot(2,3,5)
imagesc(inv(cov(data)));
title('Sample Covariance Inverse')
axis square

subplot(2,3,6)
imagesc(COVINV_EST);
title(sprintf('Regularized Covariance Inverse \\lambda=%4.2f', lambda))
axis square

% Print out the results

set(gcf,'PaperPositionMode','auto')
set(handle, 'Position', [100, 100, 1400, 900]);
print(handle, '-depsc', nameFile_EPS);
print(handle, '-djpeg', nameFile_JPG);
savefig(handle,nameFile_FIG);

end

