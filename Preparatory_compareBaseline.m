%% get the baseline firing rate for all neurons

for j = 1:length(Sum)
    fprintf('Process neuron %4.2f \n',j)
    i = j;
    event_L = sort([Sum(i).event.tsLCorr.TasteID.S; Sum(i).event.tsLCorr.TasteID.Q]);
    taste(j).Lplanning = spike2eventRasteandPSTH_NP (Sum(i).timestampN, event_L, 1000, -3000, 0);
    event_R = sort([Sum(i).event.tsRCorr.TasteID.M; Sum(i).event.tsRCorr.TasteID.O]);
    taste(j).Rplanning = spike2eventRasteandPSTH_NP (Sum(i).timestampN,event_R, 1000, -3000, 0);
   
end

%% get the index of the planning neurons
for i = 1:length(Sum)
    resp(i) = Sum(i).auROC.PlanRcorr_vsLcorr.stats.respFlag;
end
idx_p = find(resp ==1);

for i = 1:length(Sum)
    preference(i) =Sum(i).auROC.PlanRcorr_vsLcorr.p;        
end

idx_p = find(resp ==1);
idx_p_L = intersect(find(preference<0),idx_p);
idx_p_R = intersect(find(preference>0),idx_p);

%% compare the baseline firing rate with the firing during the delay

clear baseline L_firing R_firing
for i = 1:length(idx_p_L)
    L_firing(i)  = Sum(idx_p_L(i)).auROC.PlanRcorr_vsLcorr.Left.FR_avg-mean(taste(idx_p_L(i)).Lplanning.FR_avg(1:2));
    R_firing(i)  = Sum(idx_p_L(i)).auROC.PlanRcorr_vsLcorr.Right.FR_avg-mean(taste(idx_p_L(i)).Rplanning.FR_avg(1:2));
end

%%
figure;
scatter (L_firing, R_firing,'ob')
bothneg_L = intersect(find(L_firing<0),find(R_firing<0));

clear baseline L_firing_r R_firing_r
for i = 1:length(idx_p_R)
%     baseline(i)  = 1/2*mean(taste(idx_p_R(i)).Lplanning.FR_avg(1:2) + taste(idx_p_R(i)).Rplanning.FR_avg(1:2));
    L_firing_r(i)  = Sum(idx_p_R(i)).auROC.PlanRcorr_vsLcorr.Left.FR_avg-mean(taste(idx_p_R(i)).Lplanning.FR_avg(1:2));
    R_firing_r(i)  = Sum(idx_p_R(i)).auROC.PlanRcorr_vsLcorr.Right.FR_avg-mean(taste(idx_p_R(i)).Rplanning.FR_avg(1:2));
end
hold on
scatter (L_firing_r, R_firing_r,'or')
bothneg_R = intersect(find(L_firing_r<0),find(R_firing_r<0));

hold on
plot([-10,30],[-10,30],'-k')
plot([-10,30],[0,0],'--k')
plot([0,0],[-10,30],'--k')
xlabel('Left: difference in firing rates (Delay-baseline)')
ylabel('Right: difference in firing rates (Delay-baseline)')
axis square
%%
sum(abs(L_firing)<abs(R_firing))
sum(abs(L_firing_r)<abs(R_firing_r))