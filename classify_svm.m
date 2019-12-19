for j = 1:length(tastepsth)
    X =[];
    Y =[];
    for i = 1:length(taste)
        X{i} = tastepsth_bin25(j).(taste{i}).FRmatrix;
        Y{i} = i*ones([size(X{i},1),1]);
    end
    
    X = [X{1};X{2};X{3};X{4}];
    Y = [Y{1};Y{2};Y{3};Y{4}];
    
    t = templateSVM('Standardize',true,'KernelFunction','gaussian','BoxConstraint', 1);
%     t = templateSVM('Standardize',true)

    Mdl = fitcecoc(X,Y,'OptimizeHyperparameters','auto',...
        'HyperparameterOptimizationOptions',struct('AcquisitionFunctionName',...
        'expected-improvement-plus'));
    Md1 = fitcecoc(X,Y,'Learners',t, 'ClassNames',[1,2,3,4]);
    CVMdl = crossval(Mdl,'KFold',5);
    genError(j) = kfoldLoss(CVMdl);
    fprintf('Finish processing neuron # %0.f\n',j)
end