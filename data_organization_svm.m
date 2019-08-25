%% prepare data for svm decoder correct trials; data prepared from taste similarity code
%% need the tastepsth_bin25 and Errortastepsth_bin25
clear trainedClassifier confMat validationAccuracy
tasteidd = {'S', 'M', 'Q', 'O'};
testTrialN = 5;
for boot = 1:10 % rep
    for i = 1:length(tastepsth_bin25)
        for j = 1: length(tasteidd)
            spikes.(tasteidd{j}) = tastepsth_bin25(i).(tasteidd{j}).scmatrix;
            Espikes.(tasteidd{j}) = Error_tastepsth_bin25(i).(tasteidd{j}).scmatrix;
            if size(Espikes.(tasteidd{j}),1)>testTrialN
                idd = randperm(size(Espikes.(tasteidd{j}),1));
                Espikes_test.(tasteidd{j}) = Espikes.(tasteidd{j})(idd(1:testTrialN),:); % randomly choose testTrialN trials
                iddc = randperm(size(spikes.(tasteidd{j}),1));
                spikes_test.(tasteidd{j})  = spikes.(tasteidd{j})(iddc(1:testTrialN),:);
                spikes_train.(tasteidd{j}) = spikes.(tasteidd{j})(iddc(testTrialN+1:end),:);
            else
                Espikes_test.(tasteidd{j}) = Espikes.(tasteidd{j}); %
                iddc = randperm(size(spikes.(tasteidd{j}),1));
                spikes_test.(tasteidd{j})  = spikes.(tasteidd{j})(iddc(1:size(Espikes_test,1)),:);
                spikes_train.(tasteidd{j}) = spikes.(tasteidd{j})(iddc(size(Espikes_test,1)+1:end),:);
            end
        end
        
        data = [[spikes_train.S];[spikes_train.M]; [spikes_train.Q];[spikes_train.O]];
        labels = [ones(size([spikes_train.S],1),1);2*ones(size([spikes_train.M],1),1);3*ones(size([spikes_train.Q],1),1);4*ones(size([spikes_train.O],1),1)];
        [trainedClassifier{i}, validationAccuracy(i),confMat{i}] = singleUnit_svm(data,labels,1);
        
        data_test = [[spikes_test.S];[spikes_test.M]; [spikes_test.Q];[spikes_test.O]];
        label_test = [ones(size([spikes_test.S],1),1);2*ones(size([spikes_test.M],1),1);3*ones(size([spikes_test.Q],1),1);4*ones(size([spikes_test.O],1),1)];
        Edata = [[Espikes_test.S];[Espikes_test.M]; [Espikes_test.Q];[Espikes_test.O]];
        Elabels = [ones(size([Espikes_test.S],1),1);2*ones(size([Espikes_test.M],1),1);3*ones(size([Espikes_test.Q],1),1);4*ones(size([Espikes_test.O],1),1)];
        Etest(i) = length(find(predict(trainedClassifier{i},Edata)==Elabels))/length(Elabels);
        Ctest(i) = length(find(predict(trainedClassifier{i},data_test)==label_test))/length(label_test);
        
        fprintf('Finish processing neuron # %0.f\n',i)
    end
    fprintf('Finish round # %0.f\n',boot)
    Valida_acy(boot,:) = validationAccuracy;
    E_accurary(boot,:) = Etest;
    C_accurary(boot,:) = Ctest;
    clear Etest
    clear Ctest
end
%% Try to extrtact the decoding for error trials
clear trainedClassifier validationAccuracy confMat
for i = 1:length(tastepsth_bin25)
    for j = 1: length(tasteidd)
        spikes.(tasteidd{j}) = tastepsth_bin25(i).(tasteidd{j}).scmatrix;
    end
    data = [[spikes.S];[spikes.M]; [spikes.Q];[spikes.O]];
    labels = [ones(size([spikes.S],1),1);2*ones(size([spikes.M],1),1);3*ones(size([spikes.Q],1),1);4*ones(size([spikes.O],1),1)];
    [trainedClassifier{i}, validationAccuracy(i),confMat{i}] = singleUnit_svm(data,labels,1);
    fprintf('Finish processing neuron # %0.f\n',i)
end
%%
Pairline_plot([validationAccuracy',Etest'])%%
ylim([0,0.6])
%%
for i = 1:length(confMat)
    confusion(:,:,i) = confMat{i};
    temp = confusion(:,:,i);
    SM(i)            =(temp(1,2)+temp(2,1))/2;
    SO(i)            =(temp(1,4)+temp(4,1))/2;
    MO(i)            =(temp(2,4)+temp(4,2))/2;
    QO(i)            =(temp(3,4)+temp(4,3))/2;
    SQ(i)            =(temp(1,3)+temp(3,1))/2;
    MQ(i)            =(temp(2,3)+temp(3,2))/2;
end

[~,~,stats] = anova1([SM',SO',SQ',MQ',MO',QO']);
[c,~,~,gnames] = multcompare(stats);