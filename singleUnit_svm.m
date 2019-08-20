function [classificationSVM, validationAccuracy, confMat] = singleUnit_svm(data,labels,c)
% build a svm classifier for neural decoding;
%Input:
%      data: scmatrix from the psth; each row is a trial; each collum is a
%      time bin
%      labels: the id of the trial
%      parameters for controling the margin(regularization)
if nargin <3
    c = 1;
end
template = templateSVM(...
    'KernelFunction', 'linear', ...
    'PolynomialOrder', [], ...
    'KernelScale', 'auto', ...
    'BoxConstraint', c, ...
    'Standardize', true);

classificationSVM = fitcecoc(...
    data, ...
    labels, ...
    'Learners', template, ...
    'Coding', 'onevsone', ...
    'ClassNames', [1; 2; 3; 4]);

% svmPredictFcn = @(x) predict(classificationSVM, x);
% trainedClassifier.predictFcn = @(x) svmPredictFcn(predictorExtractionFcn(x));
% Add additional fields to the result struct
% trainedClassifier.ClassificationSVM = classificationSVM;
% trainedClassifier.About = 'This struct is a trained model exported from Classification Learner R2018a.';
% trainedClassifier.HowToPredict = sprintf('To make predictions on a new predictor column matrix, X, use: \n  yfit = c.predictFcn(X) \nreplacing ''c'' with the name of the variable that is this struct, e.g. ''trainedModel''. \n \nX must contain exactly 10 columns because this model was trained using 10 predictors. \nX must contain only predictor columns in exactly the same order and format as your training \ndata. Do not include the response column or any columns you did not import into the app. \n \nFor more information, see <a href="matlab:helpview(fullfile(docroot, ''stats'', ''stats.map''), ''appclassification_exportmodeltoworkspace'')">How to predict using an exported model</a>.');


% partitionedModel = crossval(trainedClassifier.ClassificationSVM, 'KFold', 10);
partitionedModel = crossval(classificationSVM, 'KFold', 10);

% Compute validation predictions
[validationPredictions, validationScores] = kfoldPredict(partitionedModel);
confMat = confusionmat(labels, validationPredictions);
confMat = confMat./repmat(sum(confMat,2),1,size(confMat,2))
% Compute validation accuracy
validationAccuracy = 1 - kfoldLoss(partitionedModel, 'LossFun', 'ClassifError');
end