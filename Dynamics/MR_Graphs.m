% This code plots the mean reversion graphs that follow from the regression
% results

clear all
close all
clc

% Input: 4 vectors that show the mean reversion over the years for TOP1,
% TOP3 and columns 2 and 3 in table 2

TOP1_C2 = [-0.016,-0.038,-0.082,-0.16,-0.34,-0.79];
TOP3_C2 = [-0.014,-0.036,-0.082,-0.16,-0.32,-0.73];
TOP1_C3 = [-0.013,-0.033,-0.069,-0.14,-0.30,-0.72];
TOP3_C3 = [-0.011,-0.031,-0.070,-0.14,-0.28,-0.67];

t = [1,2,5,10,20,50];

% Top 1 in column 2
figure(1)
plot(t,TOP1_C2,'--gs','MarkerIndices',1:6,...
    'LineWidth',2,...
    'MarkerSize',10,...
    'MarkerEdgeColor','b',...
    'MarkerFaceColor',[0.5,0.5,0.5])
xlabel('Year')
ylabel('Mean Reversion Coefficient')
title('Top 1 and Column 2')
saveas(gcf,'Results/MR_12.png')

% Top 3 in column 2
figure(2)
plot(t,TOP3_C2,'--gs','MarkerIndices',1:6,...
    'LineWidth',2,...
    'MarkerSize',10,...
    'MarkerEdgeColor','b',...
    'MarkerFaceColor',[0.5,0.5,0.5])
xlabel('Year')
ylabel('Mean Reversion Coefficient')
title('Top 3 and Column 2')
saveas(gcf,'Results/MR_32.png')

% Top 1 in column 3
figure(3)
plot(t,TOP1_C3,'--gs','MarkerIndices',1:6,...
    'LineWidth',2,...
    'MarkerSize',10,...
    'MarkerEdgeColor','b',...
    'MarkerFaceColor',[0.5,0.5,0.5])
xlabel('Year')
ylabel('Mean Reversion Coefficient')
title('Top 1 and Column 3')
saveas(gcf,'Results/MR_13.png')

% Top 3 in column 3
figure(4)
plot(t,TOP3_C3,'--gs','MarkerIndices',1:6,...
    'LineWidth',2,...
    'MarkerSize',10,...
    'MarkerEdgeColor','b',...
    'MarkerFaceColor',[0.5,0.5,0.5])
xlabel('Year')
ylabel('Mean Reversion Coefficient')
title('Top 3 and Column 3')
saveas(gcf,'Results/MR_33.png')