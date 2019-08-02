%% barplot and statistics test for barplot
function [p,h]=barplot_test(a,b);
% INPUT:
%       a is the data from lesion animal
%       b is the data from control animal
a_m=mean(a);
a_sem=std(a)/sqrt(length(a)); % standard error
b_m=mean(b);
b_sem=std(b)/sqrt(length(b)); % standard error
% name={[],'6-OHDA','Saline',[]};
% figure;
%c=bar([a_m,b_m],0.5,'EdgeColor',[0 0 0],'LineWidth',1.5);
bar(1,a_m,'EdgeColor',[0 0 0],'FaceColor','r','LineWidth',1);hold on
bar(2,b_m,'EdgeColor',[0 0 1],'FaceColor','g','LineWidth',1);
set(gca,'XTick',0:3)
% set(gca,'xticklabel',name)
hold on;
e=errorbar([a_m, b_m], [a_sem,b_sem],'LineWidth',1.5);
% plot([1 1],[mean(a_m) mean(a_m)+a_sem],'r');hold on;
e.Color ='k';
e.LineStyle = 'none';
% e.CapSize=16;
% drawnow
% errorbar_tick(e)

% c(1).FaceColor='r';
text(1,a_m/2,['n=' num2str(length(a))],'FontSize',12)
text(2,b_m/2,['n=' num2str(length(b))],'FontSize',12)
box off
set(gca,'TickDir','out')
set(gca,'fontsize',12)
set(gca,'TickLengt', [0.02 0.035]);
hold off
if kstest((a-mean(a))./std(a)) | kstest((b-mean(b))./std(b))
    [p,h]=ranksum(a,b)
    disp('not normal distribution, use rank sum test');
else
    [h,p] = ttest2(a,b);
    disp('Normal distribution, use t-test')
end
text(1.5,(b_m+a_m)/2,['p=' num2str(p)],'FontSize',12)