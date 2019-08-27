%% start to organize the file
% for delay stimulation
file = {'RVKC419\2019-05-07_RVKC419','RVKC419\2019-05-09_RVKC419','RVKC419\2019-05-11_RVKC419',...
    'RVKC430\2019-07-16_RVKC430','RVKC430\2019-07-18_RVKC430','RVKC430\2019-07-20_RVKC430',...
    'RVKC434\2019-07-25_RVKC434','RVKC434\2019-07-27_RVKC434','RVKC434\2019-07-29_RVKC434',...
    'RVKC436\2019-07-29_RVKC436','RVKC436\2019-07-31_RVKC436','RVKC436\2019-08-02_RVKC436'};
for i = 1:length(file)
    path = ['F:\Taste Discrimination\Optogenetics\TasteDiscrimination-OptoExperiment\',...
        file{i},'\'];
    data_delay(i) = performance_TAC_opto(path,1);
end
% for sampling stimulation
file = {'RVKC419\2019-05-08_RVKC419','RVKC419\2019-05-10_RVKC419','RVKC419\2019-05-12_RVKC419',...
    'RVKC430\2019-07-17_RVKC430','RVKC430\2019-07-21_RVKC430','RVKC430\2019-07-22_RVKC430',...
    'RVKC434\2019-07-24_RVKC434','RVKC434\2019-07-26_RVKC434','RVKC434\2019-07-28_RVKC434',...
    'RVKC436\2019-07-30_RVKC436','RVKC436\2019-08-01_RVKC436'};
for i = 1:length(file)
    path = ['F:\Taste Discrimination\Optogenetics\TasteDiscrimination-OptoExperiment\',...
        file{i},'\'];
    data_sampling(i) = performance_TAC_opto(path,0);
end
%%
save('summaryOptoData.mat','data_delay','data_sampling')
%% extract the sampling duration and latency for stimulation during the sampling
for i = 1:length(data_sampling)
    temp = data_sampling(i).raw_decision;
    for j = 1:size(temp,1)
        sample = temp{j,6};
        leftR  = temp{j,8};
        rightR  = temp{j,9};
        trialType = temp{j,2};
        cor_error = temp{j,3};
        sample_duration(j) = sample(end)-sample(1);
        if trialType == 1 && cor_error ==1 % left correct trial
            delay_duration(j) = leftR(1)-sample(end);
        elseif trialType == 1 && cor_error ==0 % left error trial
            delay_duration(j) = rightR(1)-sample(end);
        elseif trialType == 2 && cor_error ==1 % right correct trial
            delay_duration(j) = rightR(1)-sample(end);
        elseif trialType ==2 && cor_error ==0 % right error trial
            delay_duration(j) = leftR(1)-sample(end);
        end
    end
    data_sampling(i).sample_duration = sample_duration;
    data_sampling(i).delay_duration  = delay_duration;
    clear sample_duration delay_duration
end
%% extract the sampling duration and latency for stimulation during the delay
for i = 1:length(data_delay)
    temp = data_delay(i).raw_decision;
    for j = 1:size(temp,1)
        sample = temp{j,6};
        leftR  = temp{j,8};
        rightR  = temp{j,9};
        trialType = temp{j,2};
        cor_error = temp{j,3};
        sample_duration(j) = sample(end)-sample(1);
        if trialType == 1 && cor_error ==1 % left correct trial
            delay_duration(j) = leftR(1)-sample(end);
        elseif trialType == 1 && cor_error ==0 % left error trial
            delay_duration(j) = rightR(1)-sample(end);
        elseif trialType == 2 && cor_error ==1 % right correct trial
            delay_duration(j) = rightR(1)-sample(end);
        elseif trialType ==2 && cor_error ==0 % right error trial
            delay_duration(j) = leftR(1)-sample(end);
        end
    end
    data_delay(i).sample_duration = sample_duration;
    data_delay(i).delay_duration  = delay_duration;
    clear sample_duration delay_duration
end
%%
save('summaryOptoData.mat','data_delay','data_sampling')
%% for light stimulation sampling: parameters for both correct and error trials 
for i = 1:length(data_sampling)
    data_sampling(i).sti_avg_sample = mean(data_sampling(i).sample_duration(data_sampling(i).stiTrial));
    data_sampling(i).sti_avg_delay  = mean(data_sampling(i).delay_duration(data_sampling(i).stiTrial));
    sample = data_sampling(i).sample_duration;
    delay  = data_sampling(i).delay_duration;
    sample(data_sampling(i).stiTrial)=[];
    delay(data_sampling(i).stiTrial)=[];
    data_sampling(i).nonsti_avg_sample= mean(sample); 
    data_sampling(i).nonsti_avg_delay= mean(delay); 
end

Pairline_plot([([data_sampling.nonsti_avg_sample])', ([data_sampling.sti_avg_sample])'])
xticklabels({'','None','Light',''})
ylim([0,0.6])
ylabel('Sampling duration (s)')
title('Light Stimulation covering sampling')

Pairline_plot([([data_sampling.nonsti_avg_delay])', ([data_sampling.sti_avg_delay])'])
xticklabels({'','None','Light',''})
ylim([0,3])
ylabel('Delay (s)')
title('Light Stimulation covering sampling')
%% for light stimulation delay: parameters for both correct and error trials 
for i = 1:length(data_delay)
    data_delay(i).sti_avg_sample = mean(data_delay(i).sample_duration(data_delay(i).stiTrial));
    data_delay(i).sti_avg_delay  = mean(data_delay(i).delay_duration(data_delay(i).stiTrial));
    sample = data_delay(i).sample_duration;
    delay  = data_delay(i).delay_duration;
    sample(data_delay(i).stiTrial)=[];
    delay(data_delay(i).stiTrial)=[];
    data_delay(i).nonsti_avg_sample= mean(sample); 
    data_delay(i).nonsti_avg_delay= mean(delay); 
end

Pairline_plot([([data_delay.nonsti_avg_sample])', ([data_delay.sti_avg_sample])'])
xticklabels({'','None','Light',''})
ylim([0,0.6])
ylabel('Sampling duration (s)')
title('Light Stimulation covering delay')

Pairline_plot([([data_delay.nonsti_avg_delay])', ([data_delay.sti_avg_delay])'])
xticklabels({'','None','Light',''})
ylim([0,3])
ylabel('Delay (s)')
title('Light Stimulation covering delay')
%%
for i = 1:length(data_delay)
    sti = data_delay(i).stiTrial;
    cor_error = cell2mat(data_delay(i).raw_decision(:,3));
    cor_idx   = find(cor_error==1);
    error_idx = find(cor_error==0);
    cor_sti   = intersect(sti,cor_idx);
    error_sti = intersect(sti,error_idx);
    all_trial = 1:size(data_delay(i).raw_decision,1);
    non_sti   = setdiff(all_trial,sti);
    cor_nonsti = intersect(non_sti,cor_idx);
    error_nonsti = intersect(non_sti,error_idx);
    data_delay(i).sti_avg_sample_cor = mean(data_delay(i).sample_duration(cor_sti));
    data_delay(i).sti_avg_delay_cor  = mean(data_delay(i).delay_duration(cor_sti));
    data_delay(i).sti_avg_sample_error = mean(data_delay(i).sample_duration(error_sti));
    data_delay(i).sti_avg_delay_error  = mean(data_delay(i).delay_duration(error_sti));
    data_delay(i).nonsti_avg_sample_cor = mean(data_delay(i).sample_duration(cor_nonsti));
    data_delay(i).nonsti_avg_delay_cor  = mean(data_delay(i).delay_duration(cor_nonsti));
    data_delay(i).nonsti_avg_sample_error = mean(data_delay(i).sample_duration(error_nonsti));
    data_delay(i).nonsti_avg_delay_error  = mean(data_delay(i).delay_duration(error_nonsti));   
end

Pairline_plot([([data_delay.nonsti_avg_delay_cor])', ([data_delay.nonsti_avg_delay_error])'])
xticklabels({'','None-Correct','None-Error',''})
ylim([0,3])
ylabel('Delay (s)')
title('Light Stimulation covering delay')

Pairline_plot([([data_delay.sti_avg_delay_cor])', ([data_delay.sti_avg_delay_error])'])
xticklabels({'','Light-Correct','Light-Error',''})
ylim([0,3])
ylabel('Delay (s)')
title('Light Stimulation covering delay')
%%
save('summaryOptoData.mat','data_delay','data_sampling')
%% Let's look at the control data
%% start to organize the file
% for delay stimulation
file = {'RVKC432\2019-07-16_RVKC432','RVKC432\2019-07-18_RVKC432','RVKC432\2019-07-20_RVKC432',...
    'RVKC433\2019-07-16_RVKC433','RVKC433\2019-07-18_RVKC433','RVKC433\2019-07-20_RVKC433',...
    'RVKC437\2019-07-25_RVKC437','RVKC437\2019-07-27_RVKC437','RVKC437\2019-07-29_RVKC437',...
    'RVKC435\2019-07-29_RVKC435','RVKC435\2019-07-31_RVKC435','RVKC435\2019-08-02_RVKC435'};
for i = 1:length(file)
    path = ['F:\Taste Discrimination\Optogenetics\TasteDiscrimination-OptoExperiment\',...
        file{i},'\'];
    data_delay(i) = performance_TAC_opto(path,1);
end
% for sampling stimulation
file = {'RVKC432\2019-07-17_RVKC432','RVKC432\2019-07-19_RVKC432','RVKC432\2019-07-21_RVKC432',...
    'RVKC433\2019-07-17_RVKC433','RVKC433\2019-07-19_RVKC433','RVKC433\2019-07-21_RVKC433',...
    'RVKC437\2019-07-24_RVKC437','RVKC437\2019-07-26_RVKC437','RVKC437\2019-07-28_RVKC437',...
    'RVKC435\2019-07-30_RVKC435','RVKC435\2019-08-01_RVKC435','RVKC435\2019-08-03_RVKC435'};
for i = 1:length(file)
    path = ['F:\Taste Discrimination\Optogenetics\TasteDiscrimination-OptoExperiment\',...
        file{i},'\'];
    data_sampling(i) = performance_TAC_opto(path,0);
end
%%
for i = 1:length(data_sampling)
    temp = data_sampling(i).raw_decision;
    for j = 1:size(temp,1)
        ili(j) = mean(diff(temp{j,6}));
    end
    data_sampling(i).ili = ili;
    clear ili
end
for i = 1:length(data_sampling)
    data_sampling(i).sti_avg_ili = mean(data_sampling(i).ili(data_sampling(i).stiTrial));
    ili = data_sampling(i).ili;
    ili(data_sampling(i).stiTrial)=[];
    data_sampling(i).nonsti_avg_ili= mean(ili); 
end

Pairline_plot([([data_sampling.nonsti_avg_ili])', ([data_sampling.sti_avg_ili])'])
xticklabels({'','None','Light',''})
ylim([0,0.2])
ylabel('Inter-lick interval(s)')
title('Light Stimulation covering sampling')
