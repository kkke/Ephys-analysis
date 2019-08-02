%% barplot and statistics test for barplot
function [p,h]=barplot_test(a,b);
% INPUT:
%       a is the data from one group
%       b is the data from the other group
a_m=mean(a);
a_sem=std(a)/sqrt(length(a)); % standard error
b_m=mean(b);
b_sem=std(b)/sqrt(length(b)); % standard error
figure
hold on
bar(1,a_m,'EdgeColor',[0 0 0],'FaceColor','m','LineWidth',1);
bar(2,b_m,'EdgeColor',[0 0 1],'FaceColor','r','LineWidth',1);
set(gca,'XTick',0:3)
e=errorbar([a_m, b_m], [a_sem,b_sem],'LineWidth',1.5);
e.Color ='k';
e.LineStyle = 'none';
text(1,a_m/2,['n=' num2str(length(a))],'FontSize',12)
text(2,b_m/2,['n=' num2str(length(b))],'FontSize',12)
box off
set(gca,'TickDir','out')
set(gca,'fontsize',12)
set(gca,'TickLengt', [0.01 0.01]);
if kstest((a-mean(a))./std(a)) | kstest((b-mean(b))./std(b))
    [p,h]=ranksum(a,b)
    fprintf('Not normal distribution, use rank sum test\n');
else
    [h,p] = ttest2(a,b);
    fprintf('Normal distribution, use t-test\n')
end
text(2,(b_m)*1.2,['p=' num2str(p)],'FontSize',12)
hold off
