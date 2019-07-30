clear; close all; clc;

%% prepare dataset

load fisheriris

species_num = grp2idx(species);

% binary classification problem
X = randn(100,10);
X(:,[1,3,5,7]) = meas(1:100,:);
y = species_num(1:100);

% 80:20
rand_num = randperm(100);
X_train = X(rand_num(1:80),:);
y_train = y(rand_num(1:80),:);

X_test = X(rand_num(81:100),:);
y_test = X(rand_num(81:100),:);

%% CV partition
c = cvpartition(y_train,'k',5);

%% feature selection
opts = statset('display','iter');
fun = @(train_data, train_label, test_data, test_label)...
    sum(predict(fitcsvm(train_data, train_label, 'KernelFunction','rbf'), test_data) ~= test_label);
[fs,history] = sequentialfs(fun, X_train, y_train,'cv',c,'options',opts,'nfeatures',2);

%% Best hyperparameters

X_trian_w_best_features = X_trian(:,fs);

Md1 = fitcsvm(X_trian_w_best_features, y_train, 'KernelFunction','rbf',...
    'OptimizeHyperparameters','auto',...
    'HyperparameterOptimizationOptions',struct('AcquisitionFunctionName',...
    'expected-improvement-plus', 'ShowPlots',ture));
