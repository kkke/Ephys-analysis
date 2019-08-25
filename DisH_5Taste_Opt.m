clear
addpath('C:\Users\ke-roberto\Documents\MATLAB\Code_Ke\')
%% load raw data and extract the stimulation signals from the currect data file
thr = 2;
file = dir('*.rhd');
dataRaw = read_Intan(file.name);
[data.OptSt_on,data.OptSt_off] = Timing_onset_offset(dataRaw.analog(4,:), dataRaw.ts, thr,300,1);

%% load the summary file for the animal
cd('D:\Data_backUp\SummaryMiceDiscrimination\')
%%
load('Summary.mat')
%% Let's make sure the every analysis is based on the most recent data
 trial = summaryData(end);
 %trial = summaryData(53);

%% For planning phase
for i = 1: size(trial.raw_decision,1)
    if ~isempty(trial.raw_decision{i,6})
        central_last(i) = trial.raw_decision{i,6}(end) + trial.raw_decision{i,5};
        central_first(i) = trial.raw_decision{i,6}(1) + trial.raw_decision{i,5};
        Left =  trial.raw_decision{i,8};
        Right = trial.raw_decision{i,9};
        lateral(i) = min([Left,Right])+ trial.raw_decision{i,5};
    end
end
sampling_duration = central_last - central_first;
delay_duration    = lateral-central_last;
k=1;
%% for delay
for i = 1: size(trial.raw_decision,1)
    for j = 1:length(data.OptSt_on)
        if data.OptSt_on(j)> central_last(i)-0.2 && data.OptSt_on(j)< central_last(i)+2 % edit by Ke to add -0.2 to the central last
            sti(k) = i;
            k = k+1;
        end
    end
end
%% for detection
for i = 1: size(trial.raw_decision,1)
    for j = 1:length(data.OptSt_on)
        if data.OptSt_on(j)< central_first(i) && data.OptSt_on(j)> central_first(i)-2
            sti(k) = i;
            k = k+1;
        end
    end
end

%% Get the performance of the trials on stimulaiton
performance_sti(:,1) = cell2mat(trial.raw_decision(sti,3));
performance_sti(:,2) = cell2mat(trial.raw_decision(sti,2));
L_R = length(find(performance_sti(:,2)==1))/(length(find(performance_sti(:,2)==1))+length(find(performance_sti(:,2)==2)));
avg_per_sti = sum(performance_sti(:,1))/size(performance_sti,1);
%% get the performance of the trials without stimulation
perf = cell2mat(trial.raw_decision(:,3));
perf(sti)=[];
avg_per_non = sum(perf)/length(perf);
figure;
h1 = bar([avg_per_non,avg_per_sti],'FaceColor','w','EdgeColor','k','LineWidth',1);
ylim([0,1])
set(gca,'xticklabel',{'Non-Sti','Sti'})
set(h1,'BarWidth',0.5)

%% for unilateral stimulation
left_sti = performance_sti(find(performance_sti(:,2)==1),1);
right_sti = performance_sti(find(performance_sti(:,2)==2),1);
left_per_sti = sum(left_sti)/length(left_sti);
right_per_sti = sum(right_sti)/length(right_sti);
perf_taste = cell2mat(trial.raw_decision(:,2));
perf_taste(sti)=[];
left_non = perf(find(perf_taste==1));
right_non = perf(find(perf_taste==2));
left_per_non = sum(left_non)/length(left_non);
right_per_non = sum(right_non)/length(right_non);
figure;
plot([left_per_non, left_per_sti]);
hold on
plot([right_per_non, right_per_sti]);
legend({'left','right'})
ylim([0,1])
xlim([0,3])
ylabel('performance')



% clearvars -except delay_duration sampling_duration avg_per_non avg_per_sti
