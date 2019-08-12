function Pairline_plot(a)
figure;
plot(a','-o','Color','k')
hold on
set(gca,'XTick',0:3)
box off
set(gca,'TickDir','out')
set(gca,'fontsize',12)
set(gca,'TickLengt', [0.01 0.01]);
bar(1,mean(a(:,1)),'EdgeColor',[0 0 0],'FaceColor','m','LineWidth',1);
bar(2,mean(a(:,2)),'EdgeColor',[0 0 1],'FaceColor','r','LineWidth',1);
ylim([0.4,1])
% name = {'','No stimulation','Light stimulation',''}
% set(gca,'xticklabel',name)
[h, p] = ttest(a(:,1),a(:,2))
title('Control group')
end