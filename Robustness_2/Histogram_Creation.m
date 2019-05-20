% Create the needed histograms

load('GridOptimization_seed1_grid3')
% muT
edges = [-0.55:0.0155:1.0];
figure
hold on 
h = histogram(paramsP(:,1),edges);
line([-0.05 -0.05], ylim, 'LineWidth', 2, 'Color', 'r')
line([0.5 0.5], ylim, 'LineWidth', 2, 'Color', 'r')
title('Distribution of estimates for \mu_T')
xlabel('\mu_T')
ylabel('Frequency')
saveas(gcf,'hist_muT.png')
hold off

% sigmaT
edges = [0:0.019:1.9];
figure
hold on 
h = histogram(paramsP(:,2),edges);
line([0.5 0.5], ylim, 'LineWidth', 2, 'Color', 'r')
line([1.4 1.4], ylim, 'LineWidth', 2, 'Color', 'r')
title('Distribution of estimates for \sigma_T')
xlabel('\sigma_T')
ylabel('Frequency')
saveas(gcf,'hist_sigmaT.png')
hold off

% tau
edges = [0.6:0.026:3.2];
figure
hold on 
h = histogram(paramsP(:,3),edges);
line([1.3 1.3], ylim, 'LineWidth', 2, 'Color', 'r')
line([3.0 3.0], ylim, 'LineWidth', 2, 'Color', 'r')
title('Distribution of estimates for \tau')
xlabel('\tau')
ylabel('Frequency')
saveas(gcf,'hist_tau.png')
hold off

% theta
edges = [-0.3:0.016:1.3];
figure
hold on 
h = histogram(paramsP(:,4)*4,edges);
line([0.2 0.2], ylim, 'LineWidth', 2, 'Color', 'r')
line([0.8 0.8], ylim, 'LineWidth', 2, 'Color', 'r')
title('Distribution of estimates for \theta')
xlabel('\theta')
ylabel('Frequency')
saveas(gcf,'hist_theta.png')
hold off

% F
edges = [-0.31:0.0231:2];
figure
hold on 
h = histogram(paramsP(:,5)*4.93*0.43/5,edges);
line([0.19 0.19], ylim, 'LineWidth', 2, 'Color', 'r')
line([1.50 1.50], ylim, 'LineWidth', 2, 'Color', 'r')
title('Distribution of estimates for F')
xlabel('F')
ylabel('Frequency')
saveas(gcf,'hist_F.png')
hold off



