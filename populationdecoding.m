% add the path to the NDT so add_ndt_paths_and_init_rand_generator can be called
toolbox_basedir_name = 'E:\MATLAB\ndt.1.0.4\'
addpath(toolbox_basedir_name);
 
% add the NDT paths using add_ndt_paths_and_init_rand_generator
add_ndt_paths_and_init_rand_generator

%%
% load bp1021spk_04B_raster_data.mat

%%
% view the rasters from one neuron
% subplot(1, 2, 1)
% imagesc(~raster_data); colormap gray
% line([500 500], get(gca, 'YLim'), 'color', [1 0 0]);
% ylabel('Trials')
% xlabel('Time (ms)')
% title('rasters')
% 
% % view the PSTH for one neuron
% subplot(1, 2, 2)
% bar(sum(raster_data));
% line([500 500], get(gca, 'YLim'), 'color', [1 0 0]);
% ylabel('Number of spikes')
% xlabel('Time (ms)')
% title('PSTH')

%%
raster_file_directory_name = 'F:\Taste Discrimination\NeuralCoding\tasteData'
save_prefix_name = 'tasteData';
bin_width = 100;
step_size = 100;
 
create_binned_data_from_raster_data(raster_file_directory_name, save_prefix_name, bin_width, step_size);

%%
load tasteData_100ms_bins_100ms_sampled.mat
 
for k = 1:65
    inds_of_sites_with_at_least_k_repeats{k} = find_sites_with_k_label_repetitions(binned_labels.stimulus_ID, k);
    num_sites_with_k_repeats(k) = length(inds_of_sites_with_at_least_k_repeats{k});
end

%% 58 neurons with more than 30 repeats
% get the data with more than 30 repeats
idx = inds_of_sites_with_at_least_k_repeats{30};
binned_data = binned_data(idx);
binned_labels.stimulus_ID = binned_labels.stimulus_ID(idx);
binned_site_info.alignment_event_time = binned_site_info.alignment_event_time(idx);
save('tasteData_100ms_bins_100ms_sampled30repeats.mat','binned_data','binned_labels','binned_site_info')
% the name of the file that has the data in binned-format
%%
binned_format_file_name = 'tasteData_100ms_bins_100ms_sampled30repeats.mat';
 
% will decode the identity of which object was shown (regardless of its position)
specific_label_name_to_use = 'stimulus_ID';
 
num_cv_splits = 5;
 
ds = basic_DS(binned_format_file_name, specific_label_name_to_use, num_cv_splits)
%%
% create a feature preprocessor that z-score normalizes each neuron
 
% note that the FP objects are stored in a cell array 
% which allows multiple FP objects to be used in one analysis
%  
% the_feature_preprocessors{1} = zscore_normalize_FP;
% %%
% 
% % create the CL object
% the_classifier = max_correlation_coefficient_CL;
%%
% the_cross_validator = standard_resample_CV(ds, the_classifier, the_feature_preprocessors);
%  
% % set how many times the outer 'resample' loop is run
% % generally we use more than 2 resample runs which will give more accurate results
% % but to save time in this tutorial we are using a small number.
%  
% the_cross_validator.num_resample_runs = 10;
%%
% run the decoding analysis
% DECODING_RESULTS = the_cross_validator.run_cv_decoding;
%  
% save_file_name = 'decoding_results'
%  
% save(save_file_name, 'DECODING_RESULTS');
%%
% create the plot results object
% % note that this object takes a string in its constructor not a cell array
% plot_obj = plot_standard_results_TCT_object(save_file_name);
%  
% % put a line at the time when the stimulus was shown
% plot_obj.significant_event_times = 0;
%  
% % display the results
% plot_obj.plot_results;
%% stats
% mkdir results
% mkdir results/shuff_results/
for i = 0:5
    run_basic_decoding_shuff(i)
end
%%
result_names{1} = 'results/population_decoding_results';
plot_obj = plot_standard_results_object(result_names);
 
 
% create the names of directories that contain the shuffled data for creating null distributions
% (this is a cell array so that multiple p-values are created when comparing results)
pval_dir_name{1} = 'results/shuff_results/';
plot_obj.p_values = pval_dir_name;
 
% use data from all time bins when creating the null distribution
plot_obj.collapse_all_times_when_estimating_pvals = 1;
 
% plot the results as usual
plot_obj.plot_results;

%% 
load('population_decoding_results.mat')
performance = squeeze(mean(DECODING_RESULTS.ZERO_ONE_LOSS_RESULTS.decoding_results,2));

% mean decoding during sampling
de_sampling = mean(performance(:,16:20),2);
de_delay1   = mean(performance(:,21:30),2);
de_delay2   = mean(performance(:,31:40),2);

bar_plot_multi([de_sampling,de_delay1, de_delay2])
ylim([0,1])
ylabel('Decoding accuracy')
names = {'','Sampling', 'Delay1','Delay2',''}
xticks([0,1,2,3,4])
xticklabels(names)
[~,~,stats] = anova1([de_sampling,de_delay1, de_delay2]);
[c,~,~,gnames] = multcompare(stats);
