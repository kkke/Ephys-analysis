function control = licking_param_chemogenetic(control)
%% extract the sampling duration and latency for stimulation during the sampling
for i = 1:length(control)
    temp = control(i).raw_decision;
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
    control(i).sample_duration = sample_duration;
    control(i).delay_duration  = delay_duration;
    clear sample_duration delay_duration
end
%%
%% for light stimulation sampling: parameters for both correct and error trials 
for i = 1:length(control)
    control(i).sti_avg_sample = mean(control(i).sample_duration);
    control(i).sti_avg_delay  = mean(control(i).delay_duration);
end