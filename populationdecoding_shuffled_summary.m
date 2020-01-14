cd('E:\Taste Discrimination\NeuralCoding\allData\results_v4\results_v1')
file = dir('*.mat');
for i = 1:length(file)
    load(file(i).name);
    confusion = DECODING_RESULTS.ZERO_ONE_LOSS_RESULTS.confusion_matrix_results.confusion_matrix;
    conf_sampling = mean(confusion(:,:,16:20),3);
    conf_delay1   = mean(confusion(:,:,21:30),3);
    conf_delay2   = mean(confusion(:,:,31:40),3);
    conf_samplingSQ(i) = (conf_sampling(2,3)+ conf_sampling(3,2))/2;
    conf_delay1SQ(i)   = (conf_delay1(2,3)+conf_delay1(3,2))/2;
    conf_delay2SQ(i)   = (conf_delay1(2,3)+conf_delay1(3,2))/2;
    
    conf_samplingMO(i) = (conf_sampling(1,4)+ conf_sampling(4,1))/2;
    conf_delay1MO(i)   = (conf_delay1(1,4)+conf_delay1(4,1))/2;
    conf_delay2MO(i)   = (conf_delay1(1,4)+conf_delay1(4,1))/2;
end

shuffle_delay2Sampling_SQ = conf_delay2SQ-conf_samplingSQ;
shuffle_delay2Sampling_MO = conf_delay2MO-conf_samplingMO;
shuffle_delay2Delay1_SQ = conf_delay2SQ-conf_delay1SQ;
shuffle_delay2Delay1_MO = conf_delay2MO-conf_delay1MO;
cd('E:\Taste Discrimination\NeuralCoding\allData\results_v2')
load('population_decoding_results.mat')
confusion = DECODING_RESULTS.ZERO_ONE_LOSS_RESULTS.confusion_matrix_results.confusion_matrix;
conf_sampling = mean(confusion(:,:,16:20),3);
conf_delay1   = mean(confusion(:,:,21:30),3);
conf_delay2   = mean(confusion(:,:,31:40),3);
true_delay2Sampling_SQ = (conf_delay2(2,3)+conf_delay2(3,2))/2-(conf_sampling(2,3)+ conf_sampling(3,2))/2;
true_delay2delay1_SQ   = (conf_delay2(2,3)+conf_delay2(3,2))/2-(conf_delay1(2,3)+ conf_delay1(3,2))/2;
length(find(shuffle_delay2Sampling_SQ>true_delay2Sampling_SQ ))/1000
length(find(shuffle_delay2Delay1_SQ>true_delay2delay1_SQ ))/1000

true_delay2Sampling_MO = (conf_delay2(1,4)+conf_delay2(4,1))/2-(conf_sampling(1,4)+ conf_sampling(4,1))/2;
true_delay2delay1_MO   = (conf_delay2(1,4)+conf_delay2(4,1))/2-(conf_delay1(1,4)+ conf_delay1(4,1))/2;

length(find(shuffle_delay2Sampling_MO>true_delay2Sampling_MO ))/1000
length(find(shuffle_delay2Delay1_MO>true_delay2delay1_MO ))/1000
