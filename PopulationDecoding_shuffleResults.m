cd('E:\Taste Discrimination\NeuralCoding\CorrError_nont_prep_L\results_v3\shuff_results')
file = dir('*.mat');
for i = 1:length(file)
    load(file(i).name)
    shuffle_per(i,:) = (DECODING_RESULTS.ZERO_ONE_LOSS_RESULTS.mean_decoding_results)'
end

cd('E:\Taste Discrimination\NeuralCoding\CorrError_prep_R\results_v3\shuff_results')
file = dir('*.mat');
for i = 1:length(file)
    load(file(i).name)
    shuffle_perR(i,:) = (DECODING_RESULTS.ZERO_ONE_LOSS_RESULTS.mean_decoding_results)'
end

shuffle_per = [shuffle_per; shuffle_perR];
t = -4.45:0.1:1.45;
sem = std(shuffle_per)./sqrt(size(shuffle_per,1))
figure;
boundedline(t, mean(shuffle_per), 2.807*sem) % 99.5% confidence interval
ylim([0.4,1])
xlim([-2,1])