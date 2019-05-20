figure
hold on
plot(paramsP(:,6));
% plot(paramsP2(:,6));
title('Loss of top 100 points')
xlabel('Points')
ylabel('Loss')
legend('Grid 2', 'Grid 3')
saveas(gcf,'Loss_Comparison.png')