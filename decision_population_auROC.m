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
%% Step2-1:plot the population auROC for preparatory + taste selective response
Timecourse_decision(Sum,Idx,Idx_dir,bothnon)
title('Taste-selective preparatory response')

Timecourse_decision(Sum,Idx,Idx_dir,setdiff(1:89,bothnon)) % in total, we have 89 preparatory response
title('Exclusive preparatory response')
%% Step3: I am trying to extract the id for Right/Left Goal neuron
% clear Idx Idx_dir
for i = 1:length(Sum)
    Idx2(i) = Sum(i).auROC.Rcorr_vsLcorr.stats.respFlag;
    Idx_dir2(i) = Sum(i).auROC.Rcorr_vsLcorr.p;
end
Timecourse_decision(Sum,Idx2,Idx_dir2)
title('Goal')
%% Step4: calculate the overlap between this two group
common = intersect(find(Idx==1),find(Idx2==1));
%% Step 4: plot the correlation between planning and action
linearcorrelation(Idx_dir(find(Idx==1)), Idx_dir2(find(Idx==1)))
xlabel('Direction preference during delay')
ylabel('Diretion preference during lateral licking')
hold on
linearcorrelation(Idx_dir(find(Idx2==1)), Idx_dir2(find(Idx2==1)))
%% Step 5: plot the correlation between planning and action (both significant different)
linearcorrelation(Idx_dir(common), Idx_dir2(common))
xlim([-1.5,1.5])
ylim([-1.5,1.5])
axis square

