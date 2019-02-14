%% For the planning neurons, let try to extract the information whether they are taste selective
load('Sum2.mat')
clearvars -except Sum
%% Step 1: get the index of the planning neurons
for i = 1:length(Sum)
    resp(i) = Sum(i).auROC.PlanRcorr_vsLcorr.stats.respFlag;
end
idx_p = find(resp ==1);
%% Step2: try to get the psth of maltose  and sucrose octaacetate right trials but aligned to the lateral licks

for j = 1:length(idx_p)
    fprintf('Process neuron %4.2f \n',j)
    i = idx_p(j);
    T = Sum(i).event.tsRCorr.TasteDel;
    Tt = unique(T(:,2));
%     M = Sum(i).event.tsRCorr.TasteID.M;
%     O = Sum(i).event.tsRCorr.TasteID.O;
%     M(:,2) = 1;
%     O(:,2) = 2;
%     MO = [M;O];
%     [B, I] = sort(MO(:,1));
%     id = MO(:,2);
%     MO_idx = id(I);
    taste(j).M = find(T(:,2) ==Tt(1)); % this is just for computing, maybe M is the 2nd
    taste(j).O = find(T(:,2) ==Tt(2));
    event = Sum(i).event.tsRCorr.FLickRSpou;
    taste(j).Mplanning = spike2eventRasteandPSTH_NP (Sum(i).timestampN, event(taste(j).M), 1000, -1000, 0);
    taste(j).Oplanning = spike2eventRasteandPSTH_NP (Sum(i).timestampN, event(taste(j).O), 1000, -1000, 0);
    [taste(j).MO.p]    = Test_auROC_dISCRIMINATION  (taste(j).Mplanning.FRmatrix, taste(j).Oplanning.FRmatrix);
    
    % prepare for shuffle to compute statistic
    % put in a single vector the p values from Right and Left
    [taste(j).MO.allp] = [taste(j).Mplanning.FRmatrix;taste(j).Oplanning.FRmatrix];
    for jj =1:1000
        [taste(j).MO.shuffl.M(:,jj), idx]= datasample(taste(j).MO.allp,size(taste(j).Mplanning.FRmatrix,1),'Replace',false);
        idx1 = 1:size(taste(j).MO.allp,1);
        taste(j).MO.shuffl.O(:,jj)       =  taste(j).MO.allp(setdiff(idx1,idx));
        taste(j).MO.shuffl.p (:,jj)       = Test_auROC_dISCRIMINATION(taste(j).MO.shuffl.M(:,jj),taste(j).MO.shuffl.O(:,jj));
    end
    if taste(j).MO.p<0
        [taste(j).MO.stats.pvalue] = length(find(taste(j).MO.shuffl.p<taste(j).MO.p))/1000;
    else
        [taste(j).MO.stats.pvalue] = length(find(taste(j).MO.shuffl.p>taste(j).MO.p))/1000;
    end
    
    if taste(j).MO.stats.pvalue<0.01
        taste(j).MO.stats.respFlag=1;
    else
        taste(j).MO.stats.respFlag=0;
    end
end
%% Step3: try to get the psth of sucrose  and  quinine left trials but aligned to the lateral licks

for j = 1:length(idx_p)
    fprintf('Process neuron %4.2f \n',j)
    i = idx_p(j);
    T = Sum(i).event.tsLCorr.TasteDel;
    Tt = unique(T(:,2));
%     M = Sum(i).event.tsRCorr.TasteID.M;
%     O = Sum(i).event.tsRCorr.TasteID.O;
%     M(:,2) = 1;
%     O(:,2) = 2;
%     MO = [M;O];
%     [B, I] = sort(MO(:,1));
%     id = MO(:,2);
%     MO_idx = id(I);
    taste(j).S = find(T(:,2) ==Tt(1)); % This is just for computing, maby S is the 2nd
    taste(j).Q = find(T(:,2) ==Tt(2));
    event = Sum(i).event.tsLCorr.FLickLSpou;
    taste(j).Splanning = spike2eventRasteandPSTH_NP (Sum(i).timestampN, event(taste(j).S), 1000, -1000, 0);
    taste(j).Qplanning = spike2eventRasteandPSTH_NP (Sum(i).timestampN, event(taste(j).Q), 1000, -1000, 0);
    [taste(j).SQ.p]    = Test_auROC_dISCRIMINATION  (taste(j).Splanning.FRmatrix, taste(j).Qplanning.FRmatrix);
    
    % prepare for shuffle to compute statistic
    % put in a single vector the p values from Right and Left
    [taste(j).SQ.allp] = [taste(j).Splanning.FRmatrix;taste(j).Qplanning.FRmatrix];
    for jj =1:1000
        [taste(j).SQ.shuffl.S(:,jj), idx]= datasample(taste(j).SQ.allp,size(taste(j).Splanning.FRmatrix,1),'Replace',false);
        idx1 = 1:size(taste(j).SQ.allp,1);
        taste(j).SQ.shuffl.Q(:,jj)       =  taste(j).SQ.allp(setdiff(idx1,idx));
        taste(j).SQ.shuffl.p (:,jj)       = Test_auROC_dISCRIMINATION(taste(j).SQ.shuffl.S(:,jj),taste(j).SQ.shuffl.Q(:,jj));
    end
    if taste(j).SQ.p<0
        [taste(j).SQ.stats.pvalue] = length(find(taste(j).SQ.shuffl.p<taste(j).SQ.p))/1000;
    else
        [taste(j).SQ.stats.pvalue] = length(find(taste(j).SQ.shuffl.p>taste(j).SQ.p))/1000;
    end
    
    if taste(j).SQ.stats.pvalue<0.01
        taste(j).SQ.stats.respFlag=1;
    else
        taste(j).SQ.stats.respFlag=0;
    end
end
%% try to summarize to see which neurons show taste selectivity
for i = 1:length(taste)
    if taste(i).MO.stats.pvalue ==0 
       shu = taste(i).MO.shuffl.p;
       ok_p = length(find(taste(i).MO.p >= shu))/1000;
       oko_p= length(find(taste(i).MO.p <= shu))/1000;
       if ok_p >=0.01 & oko_p >=0.01 % if this exist, it means shuffling preference and real preference has lots in common
           fprintf('Neuron %4.2f is not good Left\n',i)
           taste(i).MO.stats.respFlag = 0;
       end
       
    elseif taste(i).SQ.stats.pvalue ==0 
       shu2 = taste(i).SQ.shuffl.p;
       ok_p = length(find(taste(i).SQ.p >= shu2))/1000;
       oko_p= length(find(taste(i).SQ.p <= shu2))/1000;
       if ok_p >=0.01 & oko_p >=0.01 % if this exist, it means shuffling preference and real preference has lots in common
           fprintf('Neuron %4.2f is not good Right\n',i)
           taste(i).SQ.stats.respFlag = 0;
       end
       
    else
    end
end
%%
clear bothnon
x=1;
for i = 1:length(taste)
       if taste(i).MO.stats.respFlag && taste(i).SQ.stats.respFlag
           fprintf('Neuron %4.2f is taste selective \n',i)
           bothnon(x) = i;
           x = x+1;
       end
end
%%
%% try to show neurons
close all
i = 3;
a = spike2eventRasteandPSTH_NP(taste(i).Mplanning.Spike,taste(i).Mplanning.Event,100,-1500,1500);
b = spike2eventRasteandPSTH_NP(taste(i).Oplanning.Spike,taste(i).Oplanning.Event,100,-1500,1500);
c = spike2eventRasteandPSTH_NP(taste(i).Splanning.Spike,taste(i).Splanning.Event,100,-1500,1500);
d = spike2eventRasteandPSTH_NP(taste(i).Qplanning.Spike,taste(i).Qplanning.Event,100,-1500,1500);

figure;
figure; plot(a.FR_avg)
hold on
plot(b.FR_avg)
plot(c.FR_avg)
plot(d.FR_avg)
legend({'M','O','S','Q'})