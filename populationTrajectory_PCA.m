%% load data
load('Sum2.mat')
%% for each neurons, calculate the psth for all tastant only correct trials were used for analysis
% binsize = 100 ms, windows: 0 to 2.5 s
data = [];
for i = 1:length(Sum)
    event.M = Sum(i).event.tsRCorr.TasteID.M;
    event.O = Sum(i).event.tsRCorr.TasteID.O;
    event.S = Sum(i).event.tsLCorr.TasteID.S;
    event.Q = Sum(i).event.tsLCorr.TasteID.Q;
    event.all = [event.M;event.O;event.S;event.O;event.Q];
    event.all = sort(event.all);
    data(i).unit_all = spike2eventRasteandPSTH_NP(Sum(i).timestampN,event.all,100,-1000,2500);
    data(i).unit_S = spike2eventRasteandPSTH_NP(Sum(i).timestampN,event.S,100,-1000,2500);
    data(i).unit_M = spike2eventRasteandPSTH_NP(Sum(i).timestampN,event.M,100,-1000,2500);
    data(i).unit_Q = spike2eventRasteandPSTH_NP(Sum(i).timestampN,event.Q,100,-1000,2500);
    data(i).unit_SO = spike2eventRasteandPSTH_NP(Sum(i).timestampN,event.O,100,-1000,2500);
    fprintf('Finish processing neuron # %0.1f\n', i )
end
%% Let's construct the firing rate matrix for each neuron
f = fieldnames(data);
for i = 1:length(Sum)
    for j = 1:length(f)
        firing_rate.(f{j})(i,:) = data(i).(f{j}).FR_avg;
    end    
end
%% PCA analysis
[coeff,score,latent,tsquared,explained,mu] = pca((firing_rate.unit_all)');
% [coeff,score,latent,tsquared,explained,mu] = pca((firing_rate.unit_all(:,11:end))');
%% apply the PCA transformation to each case
% S_traj = ((firing_rate.unit_S(:,11:end))'-mu)*coeff(:,1:3);
S_traj = (firing_rate.unit_S'-mu)*coeff(:,1:3);
Q_traj = (firing_rate.unit_Q'-mu)*coeff(:,1:3);
M_traj = (firing_rate.unit_M'-mu)*coeff(:,1:3);
SO_traj= (firing_rate.unit_SO'-mu)*coeff(:,1:3);

figure;
% plot3(S_traj(:,1),S_traj(:,2),S_traj(:,3))
% hold on
% plot3(Q_traj(:,1),Q_traj(:,2),Q_traj(:,3))
% plot3(M_traj(:,1),M_traj(:,2),M_traj(:,3))
% plot3(SO_traj(:,1),SO_traj(:,2),SO_traj(:,3))

plot(S_traj(:,1),S_traj(:,2),'-ro')
hold on
plot(Q_traj(:,1),Q_traj(:,2),'-bo')
plot(M_traj(:,1),M_traj(:,2),'-mo')
plot(SO_traj(:,1),SO_traj(:,2),'-ko')
xlabel('Principal component 1')
ylabel('principal component 2')
legend({'S','Q','M','SO'})
save('Trajectory_PCA.mat')