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
[coeff,score,latent,tsquared,explained,mu] = pca((firing_rate.unit_all(:,11:end))');
%% apply the PCA transformation to each case
S_traj = ((firing_rate.unit_S(:,11:end))'-mu)*coeff(:,1:3);
Q_traj = ((firing_rate.unit_Q(:,11:end))'-mu)*coeff(:,1:3);
M_traj = ((firing_rate.unit_M(:,11:end))'-mu)*coeff(:,1:3);
SO_traj = ((firing_rate.unit_SO(:,11:end))'-mu)*coeff(:,1:3);

S_traj = (firing_rate.unit_S'-mu)*coeff(:,1:3);
Q_traj = (firing_rate.unit_Q'-mu)*coeff(:,1:3);
M_traj = (firing_rate.unit_M'-mu)*coeff(:,1:3);
SO_traj= (firing_rate.unit_SO'-mu)*coeff(:,1:3);

figure;
plot3(S_traj(:,1),S_traj(:,2),S_traj(:,3),'-ro')
hold on
plot3(Q_traj(:,1),Q_traj(:,2),Q_traj(:,3),'-bo')
plot3(M_traj(:,1),M_traj(:,2),M_traj(:,3),'-mo')
plot3(SO_traj(:,1),SO_traj(:,2),SO_traj(:,3),'-ko')
box on
grid on
rotate3d on


plot(S_traj(:,1),S_traj(:,2),'-ro')
hold on
plot(Q_traj(:,1),Q_traj(:,2),'-bo')
plot(M_traj(:,1),M_traj(:,2),'-mo')
plot(SO_traj(:,1),SO_traj(:,2),'-ko')
xlabel('Principal component 1')
ylabel('principal component 2')
legend({'S','Q','M','SO'})
save('Trajectory_PCA.mat')
%% calculate eucleadian distance in PCA space
S_Qdist = sqrt(sum((S_traj - Q_traj).^2,2));
M_SOdist = sqrt(sum((M_traj - SO_traj).^2,2));

S_Mdist  = sqrt(sum((S_traj - M_traj).^2,2));
Q_SOdist = sqrt(sum((Q_traj - SO_traj).^2,2));

figure;
hold on
plot(data(1).unit_S.timepoint(11:end),S_Qdist,'-ro')
plot(data(1).unit_S.timepoint(11:end),M_SOdist,'-bo')
ylim([0,25])
ylabel('Eucleadian distance in PCA space')
xlabel('Time (s)')
legend({'S-Q distance','M-SO distance'})

figure;
hold on
plot(data(1).unit_S.timepoint(11:end),S_Mdist,'-mo')
plot(data(1).unit_S.timepoint(11:end),Q_SOdist,'-ko')
ylim([0,25])
ylabel('Eucleadian distance in PCA space')
xlabel('Time (s)')
legend({'S-M distance','Q-SQ distance'})
%%
for i = 1:length(data)
    data(i).SM = abs(psth_auROC_ke(data(i).unit_S.scmatrix,data(i).unit_M.scmatrix));
    data(i).SQ = abs(psth_auROC_ke(data(i).unit_S.scmatrix,data(i).unit_Q.scmatrix));
    data(i).SSO = abs(psth_auROC_ke(data(i).unit_S.scmatrix,data(i).unit_SO.scmatrix));
    data(i).MSO = abs(psth_auROC_ke(data(i).unit_M.scmatrix,data(i).unit_SO.scmatrix));
    data(i).QSO = abs(psth_auROC_ke(data(i).unit_Q.scmatrix,data(i).unit_SO.scmatrix));
    data(i).MQ = abs(psth_auROC_ke(data(i).unit_M.scmatrix,data(i).unit_Q.scmatrix));
    fprintf('Finish processing neuron # %0.1f\n', i)
end

for i = 1:length(data)
    SQ(i,:) = data(i).SQ;
    MSO(i,:) = data(i).MSO;
    SM (i,:) = data(i).SM;
    QSO(i,:) = data(i).QSO;
    SSO(i,:) = data(i).SSO;
    MQ(i,:)  = data(i).MQ;
end


similarAction = 1/2 * (SQ + MSO);
m_similarAction = mean(similarAction);
sem_similarAction = std(similarAction)/sqrt(size(similarAction,1));

D_similarAction = 1/6 * (SM + QSO + SQ + MSO + SSO + MQ);
D_similarAction = 1/4 * (SM + QSO + SSO + MQ);
D_similarAction = 1/2 * (SM + QSO);
m_DsimilarAction = mean(D_similarAction);
sem_DsimilarAction = std(D_similarAction)/sqrt(size(D_similarAction,1));

figure; 
h1 = boundedline(data(1).unit_all.timepoint, m_similarAction, sem_similarAction)
hold on
h3 = boundedline(data(1).unit_all.timepoint, m_DsimilarAction, sem_DsimilarAction,'r')
xlim([0,2.5])
ylim([0.06,0.12])
xlabel('Time (s)')
ylabel('Normalized auROC')


for i = 1:size(similarAction,2)
    [~, p(i)] = ttest2(similarAction(:,i),D_similarAction(:,i));
end
p = p*25;
idx = p<0.05;
time = data(1).unit_all.timepoint(idx);
scatter(time,0.09* ones(size(time)),'ko')


idx = find(data(1).unit_all.timepoint>=0 & data(1).unit_all.timepoint<2.5);
a_similarAct = similarAction(:,idx);
a_similarTaste = D_similarAction (:,idx);
