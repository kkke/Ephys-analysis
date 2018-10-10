%% Step1: I am trying to extract the id for Right/Left planning neuron
load('Sum2.mat')
clearvars -except Sum % load Sum2.mat (there is one neuron in Sum.mat which does not fire during delay, but captured)
for i = 1:length(Sum)
    Idx(i) = Sum(i).auROC.PlanRcorr_vsLcorr.stats.respFlag;
    Idx_dir(i) = Sum(i).auROC.PlanRcorr_vsLcorr.p;
end
%% Step2: Plot the population auROC for planning 
Timecourse_decision(Sum,Idx,Idx_dir)
title('Planning')
%% Step3: I am trying to extract the id for Right/Left Goal neuron
clear Idx Idx_dir
for i = 1:length(Sum)
    Idx(i) = Sum(i).auROC.Rcorr_vsLcorr.stats.respFlag;
    Idx_dir(i) = Sum(i).auROC.Rcorr_vsLcorr.p;
end
Timecourse_decision(Sum,Idx,Idx_dir)
title('Goal')
%%