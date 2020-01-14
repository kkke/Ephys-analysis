function run_basic_decoding_shuff_time(shuff_num)
for i = 1:shuff_num
    fprintf('Running shuffling # %0.01f\n', i)
    % the following parameters are hard coded but could be input arguments to this function
    specific_binned_labels_names = 'stimulus_ID';  % use object identity labels to decode which object was shown
    num_cv_splits = 10;     % use 30 cross-validation splits; I have 30 repeats; it is like leave one out
    % binned_data_file_name = 'tasteData_100ms_bins_100ms_sampled30repeats.mat'; % use the data that was previously binned
    binned_data_file_name = 'all_100ms_bins_100ms_sampled30repeats.mat'; % use the data that was previously binned
    cd('E:\Taste Discrimination\NeuralCoding\allData')
    load(binned_data_file_name)
    cd('E:\Taste Discrimination\NeuralCoding\allData\results_v4')
    % shuffle the time
    % totalTime = size(binned_data{1},2);
    
    shuffle = randperm(25)+15 % shuffle from 16th to 40th;
    for j = 1:length(binned_data)
        temp = binned_data{j};
        for k = 1:length(shuffle)
            binned_data{j}(:,15+k) = temp(:,shuffle(k));
            %        a(:,15+k) = binned_data{j}(:,shuffle(k));
            
        end
    end
    save('all_100ms_bins_100ms_sampled30repeats.mat','binned_data','binned_labels','binned_site_info')
    % the name of where to save the results
    save_file_name = ['results_v1/population_decoding_results',num2str(i)];
    
    
    % create the basic objects needed for decoding
    ds = basic_DS(binned_data_file_name, specific_binned_labels_names,  num_cv_splits); % create the basic datasource object
    the_feature_preprocessors{1} = zscore_normalize_FP;  % create a feature preprocess that z-score normalizes each feature
    the_classifier = max_correlation_coefficient_CL;  % select a classifier
    % the_classifier = libsvm_CL;
    the_cross_validator = standard_resample_CV(ds, the_classifier, the_feature_preprocessors);
    
    
    % if running the regular results, create the regular cross-validator
    the_cross_validator = standard_resample_CV(ds, the_classifier, the_feature_preprocessors);
    the_cross_validator.num_resample_runs = 10; % only repeat it 10 times to get the variance
    % the name of where to save the results for regular (non-shuffled) decoding results as before
    save_file_name = ['results_v1/population_decoding_results',num2str(i)];
    %      save_file_name = 'results_svm/population_decoding_results';
    
    % we will also greatly speed up the run-time of the analysis by not creating a full TCT matrix
    % (i.e., we are only training and testing the classifier on the same time bin)
    the_cross_validator.test_only_at_training_times = 1;
    
    
    % run the decoding analysis and save the results
    DECODING_RESULTS = the_cross_validator.run_cv_decoding;
    
    % save the results
    save(save_file_name, 'DECODING_RESULTS');
    
end

end