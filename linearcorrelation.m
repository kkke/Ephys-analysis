function linearcorrelation(x,y)
figure;
scatter(x,y,'ro','filled')
hold on
mdl = fitlm(x,y,'linear');
coeff = mdl.Coefficients;
rsqure = mdl.Rsquared.Ordinary;
cept = coeff.Estimate(1);
x1        = coeff.Estimate(2);
pvalue    = coeff.pValue(2);
model_x   = linspace(min(x),max(x),100);
model_y   = x1*model_x+cept;
% figure;
plot(model_x,model_y,'r')
%legend('Lesion','Control',['y = ',num2str(x1),'*x',num2str(cept)])
text(0,0.6,['R^2=',num2str(rsqure),'  p=',num2str(pvalue)])
ylim([-1,1])
xlim([-1,1])
% xlabel('Normalized FR of plan neurons')
% ylabel('Normalized FR of Goal neurons')