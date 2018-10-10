%I am trying to write a script that can generate the roc along the time.
%% Step1: I am trying to extract the id for Right/Left planning neuron
% clearvars -except Sum % load Sum2.mat (there is one neuron in Sum.mat which does not fire during delay, but captured)
% for i = 1:length(Sum)
%     Idx(i) = Sum(i).auROC.Rcorr_vsLcorr.stats.respFlag;
%     Idx_dir(i) = Sum(i).auROC.Rcorr_vsLcorr.p;
% end

function Timecourse_decision(Sum,Idx,Idx_dir)
%         Input:
%               Sum: the data saved in Sum2.mat
%               Idx: the Index of selective neuron
%               Idx: the direction of the selectivity
% For example: 
% load('Sum2.mat')
% for i = 1:length(Sum)
%     Idx(i) = Sum(i).auROC.Rcorr_vsLcorr.stats.respFlag;
%     Idx_dir(i) = Sum(i).auROC.Rcorr_vsLcorr.p;
% end
Idx_rp = find(Idx==1 & Idx_dir>0); % id for right  neuron
Idx_lp = find(Idx==1 & Idx_dir<0); % id for left  neuron
%% Step2: Calculate the right direction preference across time ((auROC-0.5)*2)
type = {'RightRew','LeftRew'};
clear PSTH_auROC_R PSTH_auROC_L
pre = -5000;
post = 2000;
for i = 1:length(Idx_rp)
    for j = 1:length(type)
        spike = Sum(Idx_rp(i)).PSTH.(type{j}).Corr.Spike;
        event = Sum(Idx_rp(i)).PSTH.(type{j}).Corr.Event;
        data.(type{j}).psth  = spike2eventRasteandPSTH_NP(spike,event,100,pre,post);
    end
    PSTH_auROC_R(i,:) = psth_auROC_ke(data.(type{2}).psth.scmatrix, data.(type{1}).psth.scmatrix);
end
m_auROC= mean(PSTH_auROC_R,1);
sem_auROC = std(PSTH_auROC_R,1)/sqrt(size(PSTH_auROC_R,1));
time = data.RightRew.psth.timepoint;
figure
h1= boundedline(time,m_auROC,sem_auROC,'r');
ylim([-0.4,0.4])
%% Step2: Calculate the Left direction preference across time ((auROC-0.5)*2)
type = {'RightRew','LeftRew'};
for i = 1:length(Idx_lp)
    for j = 1:length(type)
        spike = Sum(Idx_lp(i)).PSTH.(type{j}).Corr.Spike;
        event = Sum(Idx_lp(i)).PSTH.(type{j}).Corr.Event;
        data.(type{j}).psth  = spike2eventRasteandPSTH_NP(spike,event,100,pre,post);
    end
    PSTH_auROC_L(i,:) = psth_auROC_ke(data.(type{2}).psth.scmatrix, data.(type{1}).psth.scmatrix);
end
m_auROC= mean(PSTH_auROC_L,1);
sem_auROC = std(PSTH_auROC_L,1)/sqrt(size(PSTH_auROC_L,1));
time = data.RightRew.psth.timepoint;
hold on
h2 = boundedline(time,m_auROC,sem_auROC);
ylim([-0.3,0.3])
ylabel('Direction Planning Preference')
xlabel('Time (s)')
%% Step3:  Calculate the non-direction preference across time ((auROC-0.5)*2)
Idx_n = find(Idx ==0);
type = {'RightRew','LeftRew'};
for i = 1:length(Idx_n)
    for j = 1:length(type)
        spike = Sum(Idx_n(i)).PSTH.(type{j}).Corr.Spike;
        event = Sum(Idx_n(i)).PSTH.(type{j}).Corr.Event;
        data.(type{j}).psth  = spike2eventRasteandPSTH_NP(spike,event,100,pre,post);
    end
    PSTH_auROC_N(i,:) = psth_auROC_ke(data.(type{2}).psth.scmatrix, data.(type{1}).psth.scmatrix);
end
m_auROC= mean(PSTH_auROC_N,1);
sem_auROC = std(PSTH_auROC_N,1)/sqrt(size(PSTH_auROC_N,1));
time = data.RightRew.psth.timepoint;
hold on
h3 = boundedline(time,m_auROC,sem_auROC,'k','transparency',0.5);
ylim([-0.3,0.3])
ylabel('Direction Planning Preference')
xlabel('Time (s)')
legend([h1,h2,h3],'Right selective','Left selective','Non-selective')