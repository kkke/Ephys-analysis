function run_basic_decoding_shuff(shuff_num)
 
% the following parameters are hard coded but could be input arguments to this function
specific_binned_labels_names = 'stimulus_ID';  % use object identity labels to decode which object was shown 
num_cv_splits = 11;     % use 30 cross-validation splits; I have 30 repeats; it is like leave one out 
% binned_data_file_name = 'tasteData_100ms_bins_100ms_sampled30repeats.mat'; % use the data that was previously binned 
binned_data_file_name = 'all_100ms_bins_100ms_sampled11repeats.mat'; % use the data that was previously binned 
 
% the name of where to save the results
save_file_name = 'results_v3/population_decoding_results';
 
 
% create the basic objects needed for decoding
ds = basic_DS(binned_data_file_name, specific_binned_labels_names,  num_cv_splits); % create the basic datasource object
the_feature_preprocessors{1} = zscore_normalize_FP;  % create a feature preprocess that z-score normalizes each feature
the_classifier = max_correlation_coefficient_CL;  % select a classifier
% the_classifier = libsvm_CL;
the_cross_validator = standard_resample_CV(ds, the_classifier, the_feature_preprocessors);  
 
if shuff_num == 0
 
    'Currently running regular decoding results'
 
     % if running the regular results, create the regular cross-validator
     the_cross_validator = standard_resample_CV(ds, the_classifier, the_feature_preprocessors);  
     the_cross_validator.num_resample_runs = 50; % only repeat it 10 times to get the variance
     % the name of where to save the results for regular (non-shuffled) decoding results as before
     save_file_name = 'results_v3/population_decoding_results';
%      save_file_name = 'results_svm/population_decoding_results';

 
else
 
    'Currently running shuffled label decoding results (data for the null distribution)'
 
    ds.randomly_shuffle_labels_before_running = 1;  % randomly shuffled the labels before running
 
    % create the cross validator as before
    the_cross_validator = standard_resample_CV(ds, the_classifier, the_feature_preprocessors);  
 
    the_cross_validator.num_resample_runs = 10;  % only do 10 resample runs to save time
 
    % don't show progress to reduce visual clutter while the code is running
    the_cross_validator.display_progress.zero_one_loss = 0;  
    the_cross_validator.display_progress.resample_run_time = 0;
 
    % save the results with the appropriate name to the shuff_results/ directory
    save_file_name = ['results_v3/shuff_results/population_decoding_results_shuff_run_' num2str(shuff_num, '%03d')];
 
end
 
% we will also greatly speed up the run-time of the analysis by not creating a full TCT matrix 
% (i.e., we are only training and testing the classifier on the same time bin)
the_cross_validator.test_only_at_training_times = 1;  
 
 
% run the decoding analysis and save the results
DECODING_RESULTS = the_cross_validator.run_cv_decoding; 
 
% save the results
save(save_file_name, 'DECODING_RESULTS');
 
end