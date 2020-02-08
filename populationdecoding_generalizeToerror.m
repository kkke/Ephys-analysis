% add the path to the NDT so add_ndt_paths_and_init_rand_generator can be called
clear
toolbox_basedir_name = 'D:\MATLAB\Ephys-analysis\ndt.1.0.4\'
addpath(toolbox_basedir_name);
 
% add the NDT paths using add_ndt_paths_and_init_rand_generator
add_ndt_paths_and_init_rand_generator

%%
raster_file_directory_name = 'E:\Taste Discrimination\NeuralCoding\allData_errorTaste'

save_prefix_name = 'all';
bin_width = 100;
step_size = 100;
 
create_binned_data_from_raster_data(raster_file_directory_name, save_prefix_name, bin_width, step_size);

the_feature_preprocessors{1} = zscore_normalize_FP;  % create a feature preprocess that z-score normalizes each feature
the_classifier = max_correlation_coefficient_CL;  % select a classifier

% id_training_name = {'S','M','Q','SO'};
% id_testing_name  = {'S_r','M_r','Q_r','SO_r'};

id_string_names = {'S','M','Q','SO'};
 
for iID = 1:4   
   the_training_label_names{iID} = {id_string_names{iID}};
   the_test_label_names{iID} = {[id_string_names{iID} '_r']};
end
%%
% load nonTaste_100ms_bins_100ms_sampled.mat
%  load all_100ms_bins_100ms_sampled.mat
load all_100ms_bins_100ms_sampled30repeats.mat
binned_data_file_name = 'all_100ms_bins_100ms_sampled30repeats.mat';
num_cv_splits = 10;     % use 30 cross-validation splits; I have 30 repeats; it is like leave one out 
specific_labels_names_to_use = 'stimulus_ID';
ds = generalization_DS(binned_data_file_name, specific_labels_names_to_use,  num_cv_splits,the_training_label_names, the_test_label_names); % create the basic datasource object

the_cross_validator = standard_resample_CV(ds, the_classifier, the_feature_preprocessors);
the_cross_validator.num_resample_runs = 10; % only repeat it 10 times to get the variance
% (i.e., we are only training and testing the classifier on the same time bin)
the_cross_validator.test_only_at_training_times = 1;  
 
 
% run the decoding analysis and save the results
DECODING_RESULTS = the_cross_validator.run_cv_decoding;
%% 58 neurons with more than 30 repeats for taste responsive neurons
%  182 neurons with more than 30 repreats for all neurons
%  124 neurons with more than 30 repeats all nonTaste neurons
%  214 neurons with more than 11 repeats
% get the data with more than 30 repeats
% for error trials, with at least 10 repeats;
idx = inds_of_sites_with_at_least_k_repeats{30};
% idx = ans;
binned_data = binned_data(idx);
binned_labels.stimulus_ID = binned_labels.stimulus_ID(idx);
binned_site_info.alignment_event_time = binned_site_info.alignment_event_time(idx)
save('all_100ms_bins_100ms_sampled30repeats.mat','binned_data','binned_labels','binned_site_info','idx')

% randidx = randperm(length(idx));
% binned_data = binned_data(idx(randidx(1:23))); % matching the neurons used for taste preparatory neurons
% binned_labels.stimulus_ID = binned_labels.stimulus_ID(idx(randidx(1:23)));
% binned_site_info.alignment_event_time = binned_site_info.alignment_event_time(idx(randidx(1:23)));
% save('nonTasteData_100ms_bins_100ms_sampled40repeats.mat','binned_data','binned_labels','binned_site_info')
% save('all_100ms_bins_100ms_sampled10repeats.mat','binned_data','binned_labels','binned_site_info','idx','randidx')
% the name of the file that has the data in binned-format
%%
% binned_format_file_name = 'allData_100ms_bins_100ms_sampled30repeats.mat';
%  
% % will decode the identity of which object was shown (regardless of its position)
% specific_label_name_to_use = 'stimulus_ID';
% %  
% num_cv_splits = 5;
%  
% ds = basic_DS(binned_format_file_name, specific_label_name_to_use, num_cv_splits)

%% stats
file = {'CorrError_prep_L','CorrError_prep_R',...
    'CorrError_nont_prep_L','CorrError_nont_prep_R','CorrError_t_prep_L','CorrError_t_prep_R'};

file = {'CorrError_nont_prep_L','CorrError_nont_prep_R','CorrError_t_prep_L','CorrError_t_prep_R'};
for j = 1:length(file)
    cd(['F:\Taste Discrimination\NeuralCoding\',file{j}])
    mkdir results_v3
    mkdir results_v3/shuff_results/
    for i = 0:5
        run_basic_decoding_shuff(i)
    end
end

file = {'CorrError_nont_prep_L','CorrError_nont_prep_R','CorrError_t_prep_L','CorrError_t_prep_R'};
a =[];
for j = 1:length(file)
    cd(['F:\Taste Discrimination\NeuralCoding\',file{j},'\results_v3\'])
    load('population_decoding_results.mat')
    performance = squeeze(mean(DECODING_RESULTS.ZERO_ONE_LOSS_RESULTS.decoding_results,2));
    de_delay2   = mean(performance(:,36:45),2);
    a =[a,de_delay2];
end
b=[a(:,1);a(:,2)];
c=[a(:,3);a(:,4)];
[h,p]=ttest2(b,c);
%%
figure;
result_names{1} = 'results_v3/population_decoding_results';
plot_obj = plot_standard_results_object(result_names);
 
 
% create the names of directories that contain the shuffled data for creating null distributions
% (this is a cell array so that multiple p-values are created when comparing results)
pval_dir_name{1} = 'results_v3/shuff_results/';
plot_obj.p_values = pval_dir_name;
 
% use data from all time bins when creating the null distribution
plot_obj.collapse_all_times_when_estimating_pvals = 1;
plot_obj.p_value_alpha_level = 0.001;
% plot the results as usual
plot_obj.plot_results;
% title('taste Prep R')
%% 
load('population_decoding_results.mat')
performance = squeeze(mean(DECODING_RESULTS.ZERO_ONE_LOSS_RESULTS.decoding_results,2));

% mean decoding during sampling
de_sampling = mean(performance(:,16:20),2);
de_delay1   = mean(performance(:,21:30),2);
de_delay2   = mean(performance(:,31:40),2);


de_delay1   = mean(performance(:,36:40),2);
de_delay2   = mean(performance(:,36:45),2);

bar_plot_multi([de_sampling,de_delay1, de_delay2])
ylim([0,1])
ylabel('Decoding accuracy')
names = {'','Sampling', 'Delay1','Delay2',''}
xticks([0,1,2,3,4])
xticklabels(names)
[~,~,stats] = anova1([de_sampling,de_delay1, de_delay2]);
% [~,~,stats] = anova1([de_sampling,de_delay1, de_delay2, de_sampling2,de_delay12, de_delay22]);
figure
[c,~,~,gnames] = multcompare(stats);
%% get the confusion matrix
confusion = DECODING_RESULTS.ZERO_ONE_LOSS_RESULTS.confusion_matrix_results.confusion_matrix;
conf_sampling = mean(confusion(:,:,16:20),3);
conf_delay1   = mean(confusion(:,:,21:30),3);
conf_delay2   = mean(confusion(:,:,31:40),3);

conf_sampling = mean(confusion(:,:,21:25),3);
conf_delay1   = mean(confusion(:,:,26:35),3);
tconf_delay2   = mean(confusion(:,:,36:45),3);

figure; imagesc(conf_sampling)
caxis([15,45])
colormap(jet)
figure; imagesc(conf_delay1)
caxis([15,45])
colormap(jet)
figure; imagesc(conf_delay2./mean(sum(conf_delay2)).*100)
caxis([15,45])
colormap(jet)

figure; imagesc((conf_delay2./mean(sum(conf_delay2)).*100-tconf_delay2./mean(sum(tconf_delay2)).*100)./(conf_delay2./mean(sum(conf_delay2)).*100+tconf_delay2./mean(sum(tconf_delay2)).*100))

difConf = conf_delay2./mean(sum(conf_delay2)).*100 -tconf_delay2./mean(sum(tconf_delay2)).*100;

figure;imagesc(conf_reor_sampling)
caxis([5,75])
title('Confusion matrix during sampling')

figure;imagesc(conf_reor_delay1)
caxis([5,75])
title('Confusion matrix during delay1')

figure;imagesc(conf_reor_delay2)
caxis([5,75])
title('Confusion matrix during delay2')