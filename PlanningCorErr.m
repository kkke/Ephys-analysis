%% Step1: I am trying to extract the id for Right/Left planning neuron
load('Sum2.mat')
clearvars -except Sum % load Sum2.mat (there is one neuron in Sum.mat which does not fire during delay, but captured)
for i = 1:length(Sum)
    Idx(i) = Sum(i).auROC.PlanRcorr_vsLcorr.stats.respFlag;
    Idx_dir(i) = Sum(i).auROC.PlanRcorr_vsLcorr.p;
end
Idx_rp = find(Idx==1 & Idx_dir>0); % id for right  neuron
Idx_lp = find(Idx==1 & Idx_dir<0); % id for left  neuron
%% Trying to calculate the preference value for error trials
Idx = find(Idx ==1);
% Idx = 1:length(Sum);
%%
for i = 1:length(Idx)
    fprintf('Process neuron %4.2f \n',i)
    right = Sum(Idx(i)).event.tsRErr.FLickRSpou;
    left  = Sum(Idx(i)).event.tsLErr.FLickRSpou;
    if length(right)<10 || length(left)<10
        p(i) =NaN;
    else
        unit.right = spike2eventRasteandPSTH_NP(Sum(Idx(i)).timestampN,right,1000,-1000,0);
        unit.left  = spike2eventRasteandPSTH_NP(Sum(Idx(i)).timestampN,left,1000,-1000,0);
        p(i) = Test_auROC_dISCRIMINATION(unit.right.scmatrix,unit.left.scmatrix);
        % Using bootstrap to perform the statistic test
        com_data = [unit.right.scmatrix; unit.left.scmatrix];
        for j = 1:1000 % 1000 iteration
            [data(:,j),idx_p] = datasample(com_data,length(unit.right.scmatrix),'Replace',false);
            idx_pall          = 1:length(com_data);
            data_left(:,j)    = com_data(setdiff(idx_pall,idx_p));
            shaffle_p(j)      = Test_auROC_dISCRIMINATION(data(:,j),data_left(:,j));
        end
        if p(i)>=0 % one tail bootstrap
            stats(i).p = length(find(shaffle_p>=p(i)))/1000;
        else
            stats(i).p = length(find(shaffle_p<=p(i)))/1000;
        end
        
        if stats(i).p>0.01
            stats(i).resp = 0;
        else
            stats(i).resp = 1;
        end
        clear data data_left shaffle_p
        
    end
end
%%
iDX2 = find(~isnan(p));
errorP = p(iDX2);
corIdx = Idx(iDX2);
for i = 1:length(corIdx)
    corP(i) = Sum(corIdx(i)).auROC.PlanRcorr_vsLcorr.p;
end

figure;
scatter(corP, errorP)
hold on
% linearcorrelation(corP,errorP)
xlim([-2,2])
ylim([-2,2])
hold on
plot([0,0],[-2,2],'--k')
plot([-2,2],[0,0],'--k')

% check how many neurons show oppositive preference
% check how many neurons are taste selective
% need to load Planingtaste.mat
% need to calculate bothnon
% need to modified the following for
% clear potent
% for i = 1:length(errorP)
%     if errorP(i)*corP(i)>0
%         potent(i) = 1
%     else
%         potent(i) = 0;
%     end
% end
% idx_ddd = iDX2(find(potent ==1));
% ccc = intersect(bothnon,idx_ddd);
% clear cccc
% for i = 1:length(ccc)
%     cccc(i) = find(iDX2 == ccc(i));
% end
% hold on
% scatter(corP(cccc),errorP(cccc),'c','fill')

%% Plot the significant preference within error trials
% k = 1;
% for i = 1:length(stats)
%     if stats(i).resp ==1
%         idx_sig(k) = i;
%         k = k+1;
%     end
% end
% hold on
% iDX2 = idx_sig;
% errorP = p(iDX2);
% corIdx = Idx(iDX2);
% clear corP
% for i = 1:length(corIdx)
%     corP(i) = Sum(corIdx(i)).auROC.PlanRcorr_vsLcorr.p;
% end
% figure;
% scatter(corP, errorP,'r')
% xlim([-2,2])
% ylim([-2,2])
% hold on
% plot([0,0],[-2,2],'--k')
% plot([-2,2],[0,0],'--k')
%% try to plot taste selective neuron separately from planning neurons
% iDX2 is the index for neurons used to plot correct vs error
% bothnon is the index for neurons with significant taste selectivity
[~,taste_idx] = intersect(iDX2, bothnon);
[~,plan_idx]  = setdiff(iDX2,bothnon);
figure
scatter(corP(taste_idx), errorP(taste_idx),'b','filled')
xlim([-2,2])
ylim([-2,2])
hold on
plot([0,0],[-2,2],'--k')
plot([-2,2],[0,0],'--k')
linearcorrelation(corP(taste_idx), errorP(taste_idx))
xlabel('Corrects: planning preference')
ylabel('Errors: planning preference')
title('Taste selective neuron within direction selectiviy')

figure
scatter(corP(plan_idx), errorP(plan_idx),'r','filled')
xlim([-2,2])
ylim([-2,2])
hold on
plot([0,0],[-2,2],'--k')
plot([-2,2],[0,0],'--k')
linearcorrelation(corP(plan_idx), errorP(plan_idx))
xlabel('Corrects: planning preference')
ylabel('Errors: planning preference')
title('Neurons mainly with direction selectiviy')
%% save the data
save('ResultsForCompareCorrectErrors.mat','plan_idx','taste_idx','corP','errorP','iDX2')
%%
type = {'Left','Right'};
PSTH_auROC_R = zeros(93,70);
for i = 1:length(Idx)
    right = Sum(Idx(i)).event.tsRErr.FLickRSpou;
    left  = Sum(Idx(i)).event.tsLErr.FLickRSpou;
    if length(right)<10 || length(left)<10
        PSTH_auROC_R(i,:) =NaN;
    else
        unit.right = spike2eventRasteandPSTH_NP(Sum(Idx(i)).timestampN,right,100,-5000,2000);
        unit.left  = spike2eventRasteandPSTH_NP(Sum(Idx(i)).timestampN,left,100,-5000,2000);
        PSTH_auROC_R(i,:) = psth_auROC_ke(unit.right.scmatrix,unit.left.scmatrix);        
    end
end
%% for left trial: I mean the neural activity is left preferential
[~,left] = intersect(Idx,Idx_lp);
[~,right] = intersect(Idx,Idx_rp);
psth_left_error = PSTH_auROC_R(left,:);
psth_right_error = PSTH_auROC_R(right,:);
psth_left_error(~any(~isnan(psth_left_error),2),:)=[];
psth_right_error(~any(~isnan(psth_right_error),2),:)=[];

figure
m_auROC= mean(psth_right_error,1);
sem_auROC = std(psth_right_error,1)/sqrt(size(psth_right_error,1));
time = unit.right.timepoint;
h1= boundedline(time,m_auROC,sem_auROC,'r');
ylim([-0.4,0.4])
hold on
m_auROC= mean(psth_left_error,1);
sem_auROC = std(psth_left_error,1)/sqrt(size(psth_left_error,1));
time = unit.left.timepoint;
h2= boundedline(time,m_auROC,sem_auROC);
ylim([-0.4,0.4])
legend([h1,h2],{'Right Selective','Left Selective'})
%%
for i = 1:length(Sum)
    Idx(i) = Sum(i).auROC.PlanRcorr_vsLcorr.stats.respFlag;
    Idx_dir(i) = Sum(i).auROC.PlanRcorr_vsLcorr.p;
end
Timecourse_decision(Sum,Idx,Idx_dir)
title('Planning')