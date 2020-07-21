%%

% Extract Taste Responsivness
% pre taste 500ms
% post taste 500 ms
addpath('C:\Users\Ke Chen\Dropbox\TasteDiscriminationAnalysis_withKe')
load('Sum.mat');
%%
% Taste Responsivness calculated with Ranksum (baseline -500ms Vs. Evoked
% 500ms); firing rate value from all trials and bins (5 bins) are vectorized (see reshape)
for k = 1:size(Sum,2)
        tasteResp(k,1) = ranksum(reshape(Sum(k).PSTH.RightRew.Maltose.FRmatrix(:,5:10),1,[]),...
                                 reshape(Sum(k).PSTH.RightRew.Maltose.FRmatrix(:,11:16),1,[]));
        tasteResp(k,2) = ranksum(reshape(Sum(k).PSTH.RightRew.Octacetate.FRmatrix(:,5:10),1,[]),...
                                 reshape(Sum(k).PSTH.RightRew.Octacetate.FRmatrix(:,11:16),1,[]));
        tasteResp(k,3) = ranksum(reshape(Sum(k).PSTH.LeftRew.Sucrose.FRmatrix(:,5:10),1,[]),...
                                 reshape(Sum(k).PSTH.LeftRew.Sucrose.FRmatrix(:,11:16),1,[]));
        tasteResp(k,4) = ranksum(reshape(Sum(k).PSTH.LeftRew.Quinine.FRmatrix(:,5:10),1,[]),...
                                 reshape(Sum(k).PSTH.LeftRew.Quinine.FRmatrix(:,11:16),1,[]));
                for t=1:4
                        if tasteResp(k,t)<0.005
                            tasteFlag(k,t)=1;
                        else
                            tasteFlag(k,t)=0;
                        end
                end
end
% Taste Responsivness calculated with Ranksum (baseline -500ms Vs. Evoked
% 500ms); firing rate value from all trials and averaged across bins
for k = 1:size(Sum,2)
        tasteResp2(k,1) = ranksum(mean(Sum(k).PSTH.RightRew.Maltose.FRmatrix(:,5:10),2),...
                                  mean(Sum(k).PSTH.RightRew.Maltose.FRmatrix(:,11:16),2));
        tasteResp2(k,2) = ranksum(mean(Sum(k).PSTH.RightRew.Octacetate.FRmatrix(:,5:10),2),...
                                  mean(Sum(k).PSTH.RightRew.Octacetate.FRmatrix(:,11:16),2));
        tasteResp2(k,3) = ranksum(mean(Sum(k).PSTH.LeftRew.Sucrose.FRmatrix(:,5:10),2),...
                                  mean(Sum(k).PSTH.LeftRew.Sucrose.FRmatrix(:,11:16),2));
        tasteResp2(k,4) = ranksum(mean(Sum(k).PSTH.LeftRew.Quinine.FRmatrix(:,5:10),2),...
                                  mean(Sum(k).PSTH.LeftRew.Quinine.FRmatrix(:,11:16),2));
                for t=1:4
                        if tasteResp2(k,t)<0.01
                            tasteFlag2(k,t)=1;
                        else
                            tasteFlag2(k,t)=0;
                        end
                end
                        if tasteResp2(k,1)<0.01
                            if mean(mean(Sum(k).PSTH.RightRew.Maltose.FRmatrix(:,5:10),2))>mean(mean(Sum(k).PSTH.RightRew.Maltose.FRmatrix(:,11:16),2));
                                tasteFlag3(k,1)=1;
                            else
                                tasteFlag3(k,1)=-1;
                            end
                        end
                        if tasteResp2(k,2)<0.01
                            if mean(mean(Sum(k).PSTH.RightRew.Octacetate.FRmatrix(:,5:10),2))>mean(mean(Sum(k).PSTH.RightRew.Octacetate.FRmatrix(:,11:16),2));
                                tasteFlag3(k,2)=1;
                            else
                                tasteFlag3(k,2)=-1;
                            end
                        end
                        if tasteResp2(k,3)<0.01
                            if mean(mean(Sum(k).PSTH.LeftRew.Sucrose.FRmatrix(:,5:10),2))>mean(mean(Sum(k).PSTH.LeftRew.Sucrose.FRmatrix(:,11:16),2));
                                tasteFlag3(k,3)=1;
                            else
                                tasteFlag3(k,3)=-1;
                            end
                        end
                        if tasteResp2(k,4)<0.01
                            if mean(mean(Sum(k).PSTH.LeftRew.Quinine.FRmatrix(:,5:10),2))>mean(mean(Sum(k).PSTH.LeftRew.Quinine.FRmatrix(:,11:16),2));
                                tasteFlag3(k,4)=1;
                            else
                                tasteFlag3(k,4)=-1;
                            end
                        end
end

% 
% Taste responsivness kruscall wallis

for k = 1:size(Sum,2)
   
        Maltose(:,1)    = mean(Sum(k).PSTH.RightRew.Maltose.FRmatrix(:,11:16),2);
        Maltose(1:size(Maltose,1),2)    = 1;
        Octacet(:,1)    = mean(Sum(k).PSTH.RightRew.Octacetate.FRmatrix(:,11:16),2);
        Octacet(1:size(Octacet,1),2)    = 2;
        Sucrose(:,1)    = mean(Sum(k).PSTH.LeftRew.Sucrose.FRmatrix(:,11:16),2);
        Sucrose(1:size(Sucrose,1),2)    = 3;
        Quinine(:,1)    = mean(Sum(k).PSTH.LeftRew.Quinine.FRmatrix(:,11:16),2);
        Quinine(1:size(Quinine,1),2)    = 4;
        
        TasteMatrix   = [Maltose;Octacet;Sucrose;Quinine];
        [p(k,1),~,stats]        = kruskalwallis(TasteMatrix(:,1),TasteMatrix(:,2),'off');
        %[c,m,h,gnames] = multcompare(stats);
        if p(k,1)<0.05
              tasteFlag4(k)=1;
        else
              tasteFlag4(k)=0;            
        end
        if p(k,1)<0.01
              tasteFlag5(k)=1;
        else
              tasteFlag5(k)=0;            
        end
        clear Maltose Sucrose Quinine Octacet
end

figure(1);
y = [sum(tasteFlag)./214;sum(tasteFlag2)./214];
bar(y);legend('M','O','S','Q'); hold on;
bar(3,sum(tasteFlag4)/214);hold on;
bar(4,sum(tasteFlag5)/214);hold on;
ylim([0 .7]); ylabel('% responsive neurons');
xticklabels({'Rank sum all trials and bins','Rank sum all trials and 1 500ms bin','kw 0.05','kw 0.01'});

%% Plot PopPSTHs
%  Maltose
[idInh] = find((tasteFlag3(:,1)) < 0);
[idExc] = find((tasteFlag3(:,1)) > 0);
N =1;
for k =1:size(idExc,1)
       popM_exc(N,:) = Sum(idExc(k)).normPSTH.RightRew.Maltose;
       N=N+1;
end
clear k; clear N; N =1;
for k =1:size(idInh,1)
       popM_inh(N,:) = Sum(idInh(k)).normPSTH.RightRew.Maltose;
       N=N+1;
end
clear idInh idExc
%  Sucrose
[idInh] = find((tasteFlag3(:,3)) < 0);
[idExc] = find((tasteFlag3(:,3)) > 0);
N =1;
for k =1:size(idExc,1)
       popS_exc(N,:) = Sum(idExc(k)).normPSTH.LeftRew.Sucrose;
       N=N+1;
end
clear k; clear N; N =1;
for k =1:size(idInh,1)
       popS_inh(N,:) = Sum(idInh(k)).normPSTH.LeftRew.Sucrose;
       N=N+1;
end
clear idInh idExc
%  Octacetate
[idInh] = find((tasteFlag3(:,2)) < 0);
[idExc] = find((tasteFlag3(:,2)) > 0);
N =1;
for k =1:size(idExc,1)
       popO_exc(N,:) = Sum(idExc(k)).normPSTH.RightRew.Octacetate;
       N=N+1;
end
clear k; clear N; N =1;
for k =1:size(idInh,1)
       popO_inh(N,:) = Sum(idInh(k)).normPSTH.RightRew.Octacetate;
       N=N+1;
end
clear idInh idExc
%  Quinine
[idInh] = find((tasteFlag3(:,4)) < 0);
[idExc] = find((tasteFlag3(:,4)) > 0);
N =1;
for k =1:size(idExc,1)
       popQ_exc(N,:) = Sum(idExc(k)).normPSTH.LeftRew.Quinine;
       N=N+1;
end
clear k; clear N; N =1;
for k =1:size(idInh,1)
       popQ_inh(N,:) = Sum(idInh(k)).normPSTH.LeftRew.Quinine;
       N=N+1;
end
clear idInh idExc

figure(100) % population PSTHs for each individual taste
subplot(2,2,1);
%plot(mean(popM_exc(:,31:70))); hold on; plot(mean(popM_inh(:,31:70)));hold on;
boundedline(1:40, mean(popM_exc(:,31:70)), std(popM_exc(:,31:70))./sqrt(size(popM_exc,1)));hold on;
boundedline(1:40, mean(popM_inh(:,31:70)), std(popM_inh(:,31:70))./sqrt(size(popM_inh,1)),'r')
plot([20 20],[0.4 0.7],'--k');
legend('exc','inh','firstLick2taste');
xlabel('Time (s)'); ylabel('FRate Norm');
xticks([0 10 20 30 40]);xticklabels({'-2','-1','0','1','2'});
title('PopPSTH - Maltose');
subplot(2,2,2);
%plot(mean(popS_exc(:,31:70))); hold on; plot(mean(popS_inh(:,31:70)));hold on;
boundedline(1:40, mean(popS_exc(:,31:70)), std(popS_exc(:,31:70))./sqrt(size(popS_exc,1)));hold on;
boundedline(1:40, mean(popS_inh(:,31:70)), std(popS_inh(:,31:70))./sqrt(size(popS_inh,1)),'r')
plot([20 20],[0.4 0.7],'--k');
legend('exc','inh','firstLick2taste');
xlabel('Time (s)'); ylabel('FRate Norm');
xticks([0 10 20 30 40]);xticklabels({'-2','-1','0','1','2'});
title('PopPSTH - Sucrose');
subplot(2,2,3);
%plot(mean(popO_exc(:,31:70))); hold on; plot(mean(popO_inh(:,31:70)));hold on;
boundedline(1:40, mean(popO_exc(:,31:70)), std(popO_exc(:,31:70))./sqrt(size(popO_exc,1)));hold on;
boundedline(1:40, mean(popO_inh(:,31:70)), std(popO_inh(:,31:70))./sqrt(size(popO_inh,1)),'r')
plot([20 20],[0.4 0.7],'--k');
legend('exc','inh','firstLick2taste');
xlabel('Time (s)'); ylabel('FRate Norm');
xticks([0 10 20 30 40]);xticklabels({'-2','-1','0','1','2'});
title('PopPSTH - Octacetate');
subplot(2,2,4);
%plot(mean(popQ_exc(:,31:70))); hold on; plot(mean(popQ_inh(:,31:70)));hold on;
boundedline(1:40, mean(popQ_exc(:,31:70)), std(popQ_exc(:,31:70))./sqrt(size(popQ_exc,1)));hold on;
boundedline(1:40, mean(popQ_inh(:,31:70)), std(popQ_inh(:,31:70))./sqrt(size(popQ_inh,1)),'r')
plot([20 20],[0.4 0.7],'--k');
legend('exc','inh','firstLick2taste');
xlabel('Time (s)'); ylabel('FRate Norm');
xticks([0 10 20 30 40]);xticklabels({'-2','-1','0','1','2'});
title('PopPSTH - Quinine');

suptitle('PopPSTH - Individual Taste responsive Neurons');

figure
bigMatExc = [popM_exc;popO_exc;popS_exc;popQ_exc];
bigMatInh = [popM_inh;popO_inh;popS_inh;popQ_inh];
figure(200);
boundedline(1:40, mean(bigMatExc(:,31:70)), std(bigMatExc(:,31:70))./sqrt(size(bigMatExc,1)));hold on;
boundedline(1:40, mean(bigMatInh(:,31:70)), std(bigMatInh(:,31:70))./sqrt(size(bigMatInh,1)),'r')
plot([20 20],[0.4 0.7],'--k');
legend('exc','inh','firstLick2taste');
xlabel('Time (s)'); ylabel('FRate Norm');
xticks([0 10 20 30 40]);xticklabels({'-2','-1','0','1','2'});
title('PopPSTH - All Taste');
% popPSTH = [popM;popS;popO;popQ];
% figure(2);
% plot(mean(popPSTH));hold on;
% plot(mean(popPSTH)+std(popPSTH));hold on;
% plot(mean(popPSTH)-std(popPSTH));hold on;
%% 
%Now Look for each kruskal wallis positive neurons (tasteFlag4), to how many (numeber and quality) taste
%each neuron respond base on the rank sum
% if a neuron is flagged to have passed the kruscal wallis and at least
% respond ton one tastants, is considered taste responsive
%
H=1;
for j=1:214
     if tasteFlag4(1,j)==1 & sum(tasteFlag2(j,:))>0
         TotTasteSel(H,:)=tasteFlag2(j,:);H=H+1;
     else
     end
end

figure(3)
subplot(1,2,1);
bar(sum(TotTasteSel)/size(TotTasteSel,1),'c');
ylabel('Fraction Taste Responsive')
xticks([1 2 3 4]);xticklabels({'M','O','S','Q'});
title('Taste Responsive Tuning Taste ID');
text(1-0.35,0.05,[num2str(sum(TotTasteSel(:,1))) '/' num2str(size(TotTasteSel,1))]);
text(2-0.35,0.05,[num2str(sum(TotTasteSel(:,2))) '/' num2str(size(TotTasteSel,1))]);
text(3-0.35,0.05,[num2str(sum(TotTasteSel(:,3))) '/' num2str(size(TotTasteSel,1))]);
text(4-0.35,0.05,[num2str(sum(TotTasteSel(:,4))) '/' num2str(size(TotTasteSel,1))]);

for f = 1:size(TotTasteSel,1)
    TotTasteSel(f,5)=sum(TotTasteSel(f,1:4));
end
TotTasteSel2(1,1)=size(find(TotTasteSel(:,5)==1),1);
TotTasteSel2(1,2)=size(find(TotTasteSel(:,5)==2),1);
TotTasteSel2(1,3)=size(find(TotTasteSel(:,5)==3),1);
TotTasteSel2(1,4)=size(find(TotTasteSel(:,5)==4),1);

TotTasteSel2(2,1)=TotTasteSel2(1,1)/size(TotTasteSel,1);
TotTasteSel2(2,2)=TotTasteSel2(1,2)/size(TotTasteSel,1);
TotTasteSel2(2,3)=TotTasteSel2(1,3)/size(TotTasteSel,1);
TotTasteSel2(2,4)=TotTasteSel2(1,4)/size(TotTasteSel,1);
f = size(TotTasteSel,1);


subplot(1,2,2);
bar(TotTasteSel2(2,:),'c');
title('Taste Responsive Tuning Taste Numbers');
text(1-0.35,0.025,[num2str(TotTasteSel2(1,1)) '/' num2str(f)]);
text(2-0.35,0.025,[num2str(TotTasteSel2(1,2)) '/' num2str(f)]);
text(3-0.35,0.025,[num2str(TotTasteSel2(1,3)) '/' num2str(f)]);
text(4-0.35,0.025,[num2str(TotTasteSel2(1,4)) '/' num2str(f)]);

figure;
pie([size(TotTasteSel,1)/length(Sum),1-size(TotTasteSel,1)/length(Sum)])
ylabel('Percentage of Taste Selective Neurons')

% clear TotTasteSel;clear H;
%%
% Population decoding part
H=1;
for j=1:214
     if tasteFlag4(1,j)==1 & sum(tasteFlag2(j,:))>0
         ID(H,1)=j;H=H+1;
     else
     end
end
clear j;
%for j=1:size(ID,1)
for j=1:214
     %Pop(j).Maltose   =Sum(ID(j)).PSTH.RightRew.Maltose.FRmatrix(:,11:16);
     Pop(j).Maltose   =Sum(j).PSTH.RightRew.Maltose.FRmatrix(:,:);
     %Pop(j).Sucrose   =Sum(ID(j)).PSTH.LeftRew.Sucrose.FRmatrix(:,11:16);
     Pop(j).Sucrose   =Sum(j).PSTH.LeftRew.Sucrose.FRmatrix(:,:);
     %Pop(j).Octacetate=Sum(ID(j)).PSTH.RightRew.Octacetate.FRmatrix(:,11:16);
     Pop(j).Octacetate=Sum(j).PSTH.RightRew.Octacetate.FRmatrix(:,:);
     %Pop(j).Quinine   =Sum(ID(j)).PSTH.LeftRew.Quinine.FRmatrix(:,11:16);
     Pop(j).Quinine   =Sum(j).PSTH.LeftRew.Quinine.FRmatrix(:,:);

end
%% plot population PSTH for planning and goal
%  Planning
timepoint = (-6:0.1:2);
timepoint = timepoint +0.01/2;
timepoint(end)=[];
id =1;
id2=1;
clear idRight idLeft
for i = 1:length(Sum)
     if Sum(i).auROC.Rcorr_vsLcorr.stats.respFlag ==0
     elseif Sum(i).auROC.Rcorr_vsLcorr.p>0
         idRight(id) = i;
         id = id+1;
     else
         idLeft(id2) = i;
         id2 = id2+1;
     end     
end
figure;
% N =1;
baseline = find(timepoint>=-1 & timepoint <0);
goal     = find(timepoint>=0 & timepoint <1);
for k =1:length(idRight)
       popM_Right(k,:) = Sum(idRight(k)).normPSTH.RightRew.Corr;
       if mean(popM_Right(k,goal))> mean(popM_Right(k,baseline))
           idRight_ex(k) = 1;
       else
           idRight_ex(k) = 0;
       end
           
end
clear k; clear N; N =1;
for k =1:length(idLeft)
       popM_Left(k,:) = Sum(idLeft(k)).normPSTH.LeftRew.Corr;
       if mean(popM_Left(k,goal))> mean(popM_Left(k,baseline))
           idLeft_ex(k) = 1;
       else
           idLeft_ex(k) = 0;
       end          
end

figure(100) % population PSTHs for each individual taste
subplot(2,2,1);
%plot(mean(popM_exc(:,31:70))); hold on; plot(mean(popM_inh(:,31:70)));hold on;
idL_ex  = find(idLeft_ex ==1);
idL_inh = find(idLeft_ex ==0);
id_ex   = find(idRight_ex ==1);
id_inh  = find(idRight_ex ==0);
%%
figure
boundedline(timepoint, mean(popM_Right(id_ex,:),1), std(popM_Right(id_ex,:),1)./sqrt(size(popM_Right(id_ex,:),1)));hold on;
boundedline(timepoint, mean(popM_Right(id_inh,:),1), std(popM_Right(id_inh,:),1)./sqrt(size(popM_Right(id_inh,:),1)),'r');hold on;
rectangle('Position',[0, min(mean(popM_Right(id_inh,:),1)),1,max(mean(popM_Right(id_ex,:),1))-min(mean(popM_Right(id_inh,:),1))])
title('Right preferential neuron only')
xlim([-1,2])
%%
figure
boundedline(timepoint, mean(popM_Left(idL_ex,:),1), std(popM_Left(idL_ex,:),1)./sqrt(size(popM_Left(idL_ex,:),1)));hold on;
boundedline(timepoint, mean(popM_Left(idL_inh,:),1), std(popM_Left(idL_inh,:),1)./sqrt(size(popM_Left(idL_inh,:),1)),'r');hold on;
rectangle('Position',[0, min(mean(popM_Left(idL_inh,:),1)),1,max(mean(popM_Left(idL_ex,:),1))-min(mean(popM_Left(idL_inh,:),1))])
title('Left preferential neuron only')
xlim([-1,2])
%%
popM_exc = [];
popM_inh = [];
popM_exc =[popM_Left(idL_ex,:);popM_Right(id_ex,:)];
popM_inh =[popM_Left(idL_inh,:);popM_Right(id_inh,:)];
figure
h1 = boundedline(timepoint, mean(popM_exc,1), std(popM_exc,1)./sqrt(size(popM_exc,1)));hold on;
h2 = boundedline(timepoint, mean(popM_inh,1), std(popM_inh,1)./sqrt(size(popM_inh,1)),'r');hold on;
rectangle('Position',[0, min(mean(popM_inh,1)),1,max(mean(popM_exc,1))-min(mean(popM_inh,1))])
title('population PSTH for goal neurons')
xlim([-1,2])
%
% boundedline(1:40, mean(popM_inh(:,31:70)), std(popM_inh(:,31:70))./sqrt(size(popM_inh,1)),'r')
% plot([20 20],[0.4 0.7],'--k');
legend([h1,h2],{'exc','inh'});
xlabel('Time (s)'); ylabel('auROC');
% xticks([0 10 20 30 40]);xticklabels({'-2','-1','0','1','2'});
% title('PopPSTH - Maltose');
%% population PSTH for planning 
clearvars -except Sum
timepoint = (-6:0.1:2);
timepoint = timepoint +0.01/2;
timepoint(end)=[];
id =1;
id2=1;
clear idRight idLeft
for i = 1:length(Sum)
     if Sum(i).auROC.PlanRcorr_vsLcorr.stats.respFlag ==0
     elseif Sum(i).auROC.Rcorr_vsLcorr.p>0
         idRight(id) = i;
         id = id+1;
     else
         idLeft(id2) = i;
         id2 = id2+1;
     end     
end
figure;
% N =1;
baseline = find(timepoint>=-6 & timepoint <-5);
goal     = find(timepoint>=-1 & timepoint <0);
for k =1:length(idRight)
       popM_Right(k,:) = Sum(idRight(k)).normPSTH.RightRew.Corr;
       if mean(popM_Right(k,goal))> mean(popM_Right(k,baseline))
           idRight_ex(k) = 1;
       else
           idRight_ex(k) = 0;
       end
           
end
clear k; clear N; N =1;
for k =1:length(idLeft)
       popM_Left(k,:) = Sum(idLeft(k)).normPSTH.LeftRew.Corr;
       if mean(popM_Left(k,goal))> mean(popM_Left(k,baseline))
           idLeft_ex(k) = 1;
       else
           idLeft_ex(k) = 0;
       end          
end

figure(100) % population PSTHs for each individual taste
subplot(2,2,1);
%plot(mean(popM_exc(:,31:70))); hold on; plot(mean(popM_inh(:,31:70)));hold on;
idL_ex  = find(idLeft_ex ==1);
idL_inh = find(idLeft_ex ==0);
id_ex   = find(idRight_ex ==1);
id_inh  = find(idRight_ex ==0);
%
figure
boundedline(timepoint, mean(popM_Right(id_ex,:),1), std(popM_Right(id_ex,:),1)./sqrt(size(popM_Right(id_ex,:),1)));hold on;
boundedline(timepoint, mean(popM_Right(id_inh,:),1), std(popM_Right(id_inh,:),1)./sqrt(size(popM_Right(id_inh,:),1)),'r');hold on;
rectangle('Position',[0, min(mean(popM_Right(id_inh,:),1)),1,max(mean(popM_Right(id_ex,:),1))-min(mean(popM_Right(id_inh,:),1))])
title('Right preferential neuron only')
xlim([-4,2])
%
figure
boundedline(timepoint, mean(popM_Left(idL_ex,:),1), std(popM_Left(idL_ex,:),1)./sqrt(size(popM_Left(idL_ex,:),1)));hold on;
boundedline(timepoint, mean(popM_Left(idL_inh,:),1), std(popM_Left(idL_inh,:),1)./sqrt(size(popM_Left(idL_inh,:),1)),'r');hold on;
rectangle('Position',[0, min(mean(popM_Left(idL_inh,:),1)),1,max(mean(popM_Left(idL_ex,:),1))-min(mean(popM_Left(idL_inh,:),1))])
title('Left preferential neuron only')
xlim([-4,2])
%
popM_exc = [];
popM_inh = [];
popM_exc =[popM_Left(idL_ex,:);popM_Right(id_ex,:)];
popM_inh =[popM_Left(idL_inh,:);popM_Right(id_inh,:)];
figure
h1 = boundedline(timepoint, mean(popM_exc,1), std(popM_exc,1)./sqrt(size(popM_exc,1)));hold on;
h2 = boundedline(timepoint, mean(popM_inh,1), std(popM_inh,1)./sqrt(size(popM_inh,1)),'r');hold on;
rectangle('Position',[0, min(mean(popM_inh,1)),1,max(mean(popM_exc,1))-min(mean(popM_inh,1))])
title('population PSTH for goal neurons')
xlim([-4,2])
%
% boundedline(1:40, mean(popM_inh(:,31:70)), std(popM_inh(:,31:70))./sqrt(size(popM_inh,1)),'r')
% plot([20 20],[0.4 0.7],'--k');
legend([h1,h2],{'exc','inh'});
xlabel('Time (s)'); ylabel('auROC');
%%
clearvars -except Sum
timepoint = (-6:0.1:2);
timepoint = timepoint +0.01/2;
timepoint(end)=[];
for i = 1:length(Sum)
    popM_left(i,:) = Sum(i).normPSTH.LeftRew.Corr;
    popM_right(i,:) = Sum(i).normPSTH.RightRew.Corr;
end
figure
h1 = boundedline(timepoint, mean(popM_left,1), std(popM_left,1)./sqrt(size(popM_left,1)));hold on;
figure
h2 = boundedline(timepoint, mean(popM_right,1), std(popM_right,1)./sqrt(size(popM_right,1)),'r');hold on;
rectangle('Position',[0, min(mean(popM_inh,1)),1,max(mean(popM_exc,1))-min(mean(popM_inh,1))])
% % Decoder part

% Neur   = 2;
% trNumM = size(Sum(Neur).PSTH.RightRew.Maltose.spikeraster,2);
% trNumO = size(Sum(Neur).PSTH.RightRew.Octacetate.spikeraster,2);
% trNumS = size(Sum(Neur).PSTH.LeftRew.Sucrose.spikeraster,2);
% trNumQ = size(Sum(Neur).PSTH.LeftRew.Quinine.spikeraster,2);
% 
% trNum = min([trNumO;trNumM;trNumS;trNumQ]);
%     
% 
%  % extract the spike raster for each tastants
%     for t =1:trNum
%         spike.S(t).spk     = Sum(Neur).PSTH.RightRew.Maltose.spikeraster(t).times;
%         spike.Q(t).spk     = Sum(Neur).PSTH.LeftRew.Quinine.spikeraster(t).times;
%         spike.M(t).spk     = Sum(Neur).PSTH.RightRew.Maltose.spikeraster(t).times;
%         spike.O(t).spk     = Sum(Neur).PSTH.RightRew.Octacetate.spikeraster(t).times;
%         %spike.W(t).spk     = psth.W.spikeraster(t).times;
%     end
%     
%     spike.S=spike.S';
%     spike.Q=spike.Q';
%     spike.M=spike.M';
%     spike.O=spike.O';
%     %spike.W=spike.W';
% 
%    spikes=spike;
% 
%  % 'spike' has fields: 'S', 'D', 'Q', 'C'-> trig
%     trig=fieldnames(spike);
%     nev=numel(trig);
%     bin = 100;
%     % decoding parameters
%     kind='units';
%     xval=1; % number of leave-out trials for each x-validation run
%     binsize=bin/1000; % size of decoding window
%     window=[0 0.5]; % total trial interval to be decoded aligned at t=0.
%     nboot=1000; % number of bootstrap runs for each shuffled stimulus
%     filesave_dec='decoding'; % prefix of files to be saved
%     
%     
%         %fun_decode_units_selectivity;
%     % OTHER PARAMETERS FOR DECODING
%     nbag=1; % for bagging, recommended if the number of trials is < 10
%     shuffled=1; % 1 to include shuffle analysis
%     fieldNames={'binsize','window','xval','nboot','nbag','shuffled','fieldNames'};
%     options=v2struct(fieldNames);
%     % initialize empty variables
%     results.all=struct('data',[],...
%         'nrun',[],'ntrain',[],'ntest',xval,'nbag',nbag,'xbins',[],'shuffled',[],'pvalue',[],'pvalue_stimulus',[]);
%     
%     % find smallest number of trials across events
%     ntrials=zeros(1,numel(trig));
%     for E=1:numel(trig)
%         [ntrials(E), nunits]=size(spikes.(trig{E}));
%     end
%     options.mintrials=min(ntrials); % pick smallest number of trials across conditions for decoding\
%     
%     tic
%     spikes_train=spikes; spikes_test=spikes;
%     results.all=decoder_multi_bag_x_selectivity(spikes_train,spikes_test,options);
%     toc
%     numfig=2;
%     fun_decode_units_selectivity_plot;