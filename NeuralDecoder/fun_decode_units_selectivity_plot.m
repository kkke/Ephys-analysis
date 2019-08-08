% plot
colors=distinguishable_colors(10);
fileload=[filesave_dec '_selectivity.mat'];
datalat=load(fileload,'results');
nstim=nev; % number of stimuli (needed below for multiple comparison across tastes)
% plot
pvalue=0.05/nstim;
figure(numfig); numfig=numfig+1; clf; hold on;
line([0 nstim],[100/nstim, 100/nstim],'linestyle','--');
accuracy=NaN(1,nstim);
accuracy(:)=100*diag(datalat.results.all.data.confusion);
for st=1:nstim
    her(st)=errorbar(st,accuracy(1,st),0,0);
    set(her(st),'Color',colors(st,:),'Marker','o','MarkerSize',10,...
        'MarkerEdgeColor',[.7 .7 .7],'MarkerFaceColor' ,colors(st,:) );
end
fprintf('Selectivity by decoding: ');
for st=1:nstim; fprintf('%s=%0.03g%%, ',trig{st},accuracy(1,st)); end
fprintf('\n                pvalues: ');
pval=datalat.results.all.pvalue_stimulus;
for st=1:nstim; fprintf('p(%s)=%0.03g, ',trig{st},pval(st)); end
indsig=find(pval<0.05/nstim);
if ~isempty(indsig)
    fprintf('\n---> %d tastes decoded (p<0.05/nstim): ',numel(indsig));
    for st=indsig; fprintf('%s, ',trig{st}); end
else
    fprintf('\n---> Neurons is not selective');
end
fprintf('\n');
% shuffled
accuracyshuff=100*diag(datalat.results.all.shuffled.confusion);
ci=100*datalat.results.all.shuffled.accuracy_ci_perstimulus;
her(nstim+1)=errorbar(1:nstim,accuracyshuff,accuracyshuff-ci(1:nstim,1),ci(1:nstim,2)-accuracyshuff);
set(her(nstim+1),'Color','k','Marker','o','MarkerSize',5,...
    'MarkerEdgeColor','k','MarkerFaceColor' , [.7 .7 .7]);
legend(her,[trig','shuffled']);
ylim([0,100]); xlim([0.5 nstim+0.5]);
set(gca,'xtick',1:numel(trig),'xticklabels',trig);
%figset(gca,'','Decoding accuracy (%)','',15);
%filename=[filesave_dec 'UnitDecodingAverage.pdf'];
%saveas(gcf,filename,'pdf');
%fprintf('Comparison plot saved to %s\n',filename);
