cd('E:\Taste Discrimination\NeuralCoding\allData_CorrTaste_downsampled_R_realign\results_v3')
load('population_decoding_results.mat')
performance_rc = squeeze(mean(DECODING_RESULTS.ZERO_ONE_LOSS_RESULTS.decoding_results,2));
cd('E:\Taste Discrimination\NeuralCoding\allData_CorrTaste_downsampled_L_realign\results_v3')
load('population_decoding_results.mat')
performance_lc = squeeze(mean(DECODING_RESULTS.ZERO_ONE_LOSS_RESULTS.decoding_results,2));

performance_c = [performance_rc;performance_lc];


cd('E:\Taste Discrimination\NeuralCoding\allData_erroTaste_R_realign\results_v3')
load('population_decoding_results.mat')
performance_re = squeeze(mean(DECODING_RESULTS.ZERO_ONE_LOSS_RESULTS.decoding_results,2));
cd('E:\Taste Discrimination\NeuralCoding\allData_erroTaste_L_realign\results_v3')
load('population_decoding_results.mat')
performance_le = squeeze(mean(DECODING_RESULTS.ZERO_ONE_LOSS_RESULTS.decoding_results,2));

performance_e = [performance_re;performance_le];
time = -1.95:0.1:0.95;
figure; 
h1 = boundedline(time,mean(performance_c), 2.8*std(performance_c)./sqrt(20),'alpha')
hold on;
h2 = boundedline(time,mean(performance_e), 2.8*std(performance_e)./sqrt(20),'r','alpha')
legend([h1,h2],{'Correct','Error'})
ylim([0.1,0.9])
ylabel('Classification accuracy')
xlabel('Time')