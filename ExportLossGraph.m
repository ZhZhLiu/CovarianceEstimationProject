function [] = ExportLossGraph(data, SIG_TRUE, lambdas, eLosses, qLosses,rootPath)

nLam = length(lambdas);
fileName = 'LossVariousLambdas';

% Handling the directory and folder

epsPath = strcat(rootPath, 'LossGraph/EPS/');
figPath = strcat(rootPath, 'LossGraph/FIG/');
jpgPath = strcat(rootPath, 'LossGraph/JPG/');

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

SAMPLE_COV = cov(data);
[bcEntropy, bcQuadratic] = getLoss(SAMPLE_COV, SIG_TRUE);

handle = figure(5);
subplot(1,2,1);
hold on
plot(lambdas,ones(nLam,1) * bcEntropy, ':r');
plot(lambdas,eLosses);
hold off
xlabel('\lambda');
title('Entropy Losses');

subplot(1,2,2);
title('Quadratic Losses');
hold on
plot(lambdas,ones(nLam,1) * bcQuadratic, ':r');
plot(lambdas,qLosses);
hold off
xlabel('\lambda');

% Export the graph
set(gcf,'PaperPositionMode','auto')
set(handle, 'Position', [100, 100, 1400, 450]);
print(handle, '-depsc', nameFile_EPS);
print(handle, '-djpeg', nameFile_JPG);
savefig(handle,nameFile_FIG);

end

