%%
clear
load('data_prep.mat');
sumData = [Data_t_prep, Data_nont_prep];
sumData = Sum_1;
clearvars -except Sum Sum_1,
%%
for i = 1:length(sumData)
    idx = find(sumData(1).PSTH.RightRew.Corr.timepoint>-1 & sumData(1).PSTH.RightRew.Corr.timepoint<0);
    RCmeanFire(i) = mean(sumData(i).PSTH.RightRew.Corr.FR_avg(idx));
    if length(sumData(i).PSTH.RightRew.Err.FR_avg)<2
        REmeanFire(i) =NaN
    else
        REmeanFire(i) = mean(sumData(i).PSTH.RightRew.Err.FR_avg(idx));
    end
    
    RCmaxFire(i) = max(sumData(i).PSTH.RightRew.Corr.FR_avg(idx));
    if length(sumData(i).PSTH.RightRew.Err.FR_avg)<2
       REmaxFire(i) =NaN
    else
       REmaxFire(i) = max(sumData(i).PSTH.RightRew.Err.FR_avg(idx));
    end
    
    LCmeanFire(i) = mean(sumData(i).PSTH.LeftRew.Corr.FR_avg(idx));
    if length(sumData(i).PSTH.LeftRew.Err.FR_avg)<2
        LEmeanFire(i) =NaN
    else
        LEmeanFire(i) = mean(sumData(i).PSTH.LeftRew.Err.FR_avg(idx));
    end
    
     LCmaxFire(i) = max(sumData(i).PSTH.LeftRew.Corr.FR_avg(idx));
    if length(sumData(i).PSTH.LeftRew.Err.FR_avg)<2
       LEmaxFire(i) =NaN
    else
       LEmaxFire(i) = max(sumData(i).PSTH.LeftRew.Err.FR_avg(idx));
    end   
    
end
%%
idx = find(isnan(REmeanFire));
idx2 = find(isnan(LEmeanFire));
RCmeanFire(idx)=[];
REmeanFire(idx)=[];
RCmaxFire(idx)=[];
REmaxFire(idx)=[];
LCmeanFire(idx)=[];
LEmeanFire(idx)=[];
LCmaxFire(idx)=[];
LEmaxFire(idx)=[];
linearcorrelation(RCmeanFire,REmeanFire)
xlim([0,20])
ylim([0,20])
linearcorrelation(RCmaxFire,REmaxFire)
xlim([0,20])
ylim([0,20])

linearcorrelation(LCmeanFire,LEmeanFire)
xlim([0,20])
ylim([0,20])
linearcorrelation(LCmaxFire,LEmaxFire)
xlim([0,20])
ylim([0,20])

linearcorrelation(LCmaxFire,REmaxFire)
xlim([0,20])
ylim([0,20])

linearcorrelation(RCmeanFire,LEmeanFire)
xlim([0,20])
ylim([0,20])
%%
figure;
scatter((RCmeanFire-REmeanFire)./(RCmeanFire+REmeanFire),(LCmeanFire-LEmeanFire)./(LCmeanFire+LEmeanFire))
linearcorrelation((RCmeanFire-REmeanFire)./(RCmeanFire+REmeanFire),(LCmeanFire-LEmeanFire)./(LCmeanFire+LEmeanFire))