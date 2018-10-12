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
    right = Sum(Idx(i)).event.tsRErr.FLickRSpou;
    left  = Sum(Idx(i)).event.tsLErr.FLickRSpou;
    if length(right)<10 || length(left)<10
        p(i) =NaN;
    else
        unit.right = spike2eventRasteandPSTH_NP(Sum(Idx(i)).timestampN,right,1000,-1000,0);
        unit.left  = spike2eventRasteandPSTH_NP(Sum(Idx(i)).timestampN,left,1000,-1000,0);
        p(i) = Test_auROC_dISCRIMINATION(unit.right.scmatrix,unit.left.scmatrix);
        % Using bootstrap to perform the statistic test
%         com_data = [unit.right.scmatrix; unit.left.scmatrix];
%         for j = 1:1000 % 1000 iteration
%             [data(:,j),idx_p] = datasample(com_data,length(unit.right.scmatrix),'Replace',false);
%             idx_pall          = 1:length(com_data);
%             data_left(:,j)    = com_data(setdiff(idx_pall,idx_p));
%             shaffle_p(j)      = Test_auROC_dISCRIMINATION(data(:,j),data_left(:,j));
%         end
%         if p(i)>=0 % one tail bootstrap
%             stats(i).p = length(find(shaffle_p>=p(i)))/1000;
%         else
%             stats(i).p = length(find(shaffle_p<=p(i)))/1000;
%         end
%         
%         if stats(i).p>0.01
%             stats(i).resp = 0;
%         else
%             stats(i).resp = 1;
%         end
%         clear data data_left shaffle_p
        
    end
end
iDX2 = find(~isnan(p));
errorP = p(iDX2);
corIdx = Idx(iDX2);
for i = 1:length(corIdx)
    corP(i) = Sum(corIdx(i)).auROC.PlanRcorr_vsLcorr.p;
end

figure;
scatter(corP, errorP)
hold on
linearcorrelation(corP,errorP)
xlim([-2,2])
ylim([-2,2])
hold on
plot([0,0],[-2,2],'--k')
plot([-2,2],[0,0],'--k')

%% 
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
% % figure;
% scatter(corP, errorP,'r')
% xlim([-2,2])
% ylim([-2,2])
% hold on
% plot([0,0],[-2,2],'--k')
% plot([-2,2],[0,0],'--k')

