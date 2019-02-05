% Create the needed histograms

load('GridOptimization_seed1_grid1')
% muT
edges = [-1.0:0.01:1.0];
figure
hold on 
h = histogram(paramsP(:,1),edges);
line([-0.5 -0.5], ylim, 'LineWidth', 2, 'Color', 'r')
line([0.5 0.5], ylim, 'LineWidth', 2, 'Color', 'r')
title('Distribution of estimates for \mu_T')
xlabel('\mu_T')
ylabel('Frequency')
saveas(gcf,'hist_muT.png')
hold off

% sigmaT
edges = [0.5:0.008:1.3]
figure
hold on 
h = histogram(paramsP(:,2),edges);
line([1.0 1.0], ylim, 'LineWidth', 2, 'Color', 'r')
line([1.8 1.8], ylim, 'LineWidth', 2, 'Color', 'r')
title('Distribution of estimates for \sigma_T')
xlabel('\sigma_T')
ylabel('Frequency')
saveas(gcf,'hist_sigmaT.png')
hold off

% tau
edges = [0.6:0.011:2.7]
figure
hold on 
h = histogram(paramsP(:,3),edges);
line([1.1 1.1], ylim, 'LineWidth', 2, 'Color', 'r')
line([2.2 2.2], ylim, 'LineWidth', 2, 'Color', 'r')
title('Distribution of estimates for \tau')
xlabel('\tau')
ylabel('Frequency')
saveas(gcf,'hist_tau.png')
hold off

% theta
edges = [3.3:0.01:4.8]
figure
hold on 
h = histogram(paramsP(:,4)*4,edges);
line([3.8 3.8], ylim, 'LineWidth', 2, 'Color', 'r')
line([4.8 4.8], ylim, 'LineWidth', 2, 'Color', 'r')
title('Distribution of estimates for \theta')
xlabel('\theta')
ylabel('Frequency')
saveas(gcf,'hist_theta.png')
hold off

% F
edges = [-0.2880:0.01:3.0439]
figure
hold on 
h = histogram(paramsP(:,5)*4.93*0.43/5,edges);
line([0.212 0.212], ylim, 'LineWidth', 2, 'Color', 'r')
line([2.5439 2.5439], ylim, 'LineWidth', 2, 'Color', 'r')
title('Distribution of estimates for F')
xlabel('F')
ylabel('Frequency')
saveas(gcf,'hist_F.png')
hold off



