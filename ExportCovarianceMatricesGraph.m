function [] = ExportCovarianceMatricesGraph(SIG_TRUE, COV_EST, SIGINV, COVINV_EST, data, fileName, lambda)

handle = figure(1);

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

nameFile_EPS = sprintf('EstimationFigures/Eps/%s.eps', fileName);
nameFile_FIG = sprintf('EstimationFigures/Fig/%s.fig', fileName);
nameFile_JPG = sprintf('EstimationFigures/Jpg/%s.jpg', fileName);

set(gcf,'PaperPositionMode','auto')
set(handle, 'Position', [100, 100, 1400, 900]);
print(handle, '-depsc', nameFile_EPS);
print(handle, '-djpeg', nameFile_JPG);
savefig(handle,nameFile_FIG);

end

