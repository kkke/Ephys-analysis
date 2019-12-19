%% extract stimuli
function sumData = performance_TAC_opto(path,delay_idx)
% This function is modified by Ke to be compatible with 4 tastant
% discrimination
cd(path)
load('data.mat');
tasteID =[data.tastant_1.line,data.tastant_2.line,data.tastant_3.line,data.tastant_4.line];
x       = {data.tastant_1.id,data.tastant_2.id,data.tastant_3.id,data.tastant_4.id};
cwd = pwd;

LeftTaste=[];
for i=1:length(data.leftID)
    if ~iscell(data.leftID)
        LeftTaste=data.(data.leftID).dig;
        LeftTaste(:,2)=1;
        LeftTaste(:,3)=data.(data.leftID).line;
    else
        Left=data.(data.leftID{i}).dig;
        Left(:,2)=1; % use 1 for taste associated with left spout
        Left(:,3)=data.(data.leftID{i}).line;
        LeftTaste=[LeftTaste; Left];
        clear Left
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
RightTaste=[];
for i=1:length(data.rightID)
    if ~iscell(data.rightID)
        RightTaste=data.(data.rightID).dig;
        RightTaste(:,2)=2; % use 3 for taste associated with right spout
        RightTaste(:,3)=data.(data.rightID).line;
    else
        Right=data.(data.rightID{i}).dig;
        Right(:,2)=2; % use 3 for taste associated with right spout
        Right(:,3)=data.(data.rightID{i}).line;
        RightTaste=[RightTaste; Right];
        clear Right
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
stimuli         =[LeftTaste;RightTaste];
%% sort stimuli
[~, Stimuli_ind]=sort(stimuli(:,1)); % sort the event (W and NaCl) as ascending
for i=1:length (Stimuli_ind)
    st(i,:)=stimuli(Stimuli_ind(i),:);
end
clear stimuli
stimuli=st;
%%
%sumData.id=filename;
sumData.mouseID = data.mouseID;
sumData.date    = data.date   ;
sumData.stimuli=stimuli;
%%

% LeftCorrect=data.Left_Correct.dig;
% RightCorrect=data.Right_Correct.dig;
decision_onset=data.Up.dig_offset;

%% 1st, find out the trial with no licks to the centeral port
trial_noSampling=[];
% for i=1:size(st,1)
%     if isempty(find(decision_onset-st(i,1)>2 & decision_onset-st(i,1)<7 ))
%         trial_noSampling=[trial_noSampling i];
%     end
% end
for i=1:size(st,1)
    if isempty(find(decision_onset-st(i,1)>0 & decision_onset-st(i,1)<5 ))
        trial_noSampling=[trial_noSampling i];
    end
end
sumData.noSampling=trial_noSampling;
%%
st(trial_noSampling,:)=[];

if size(st,1)~= length(decision_onset)
    warning('Be careful to check the data, Probably there is something wrong')
end
trial_end=data.Down.dig;
if size(st,1)~= length(trial_end)
    warning('Be careful to check the data, Probably there is something wrong')
    if size(st,1)<length(trial_end)% it means recorded more Up/Down lateral port signal
        for i=1:size(st,1)
            ind=find(trial_end>st(i,1));
            trial_end(i)=trial_end(ind(1));
        end
    else
        st(end,:)=[];   % it means recorded less Up/Down lateral port signal (likely happended for the last trial)
    end
end
trial_noResponse=[];
LeftError.trial=[]
RightError.trial=[];
for i=1:size(st,1)
    raw{i,1}=st(i,2);
    raw{i,2}=st(i,1);
    c=spike2eventRasteandPSTH_NP(data.centSp.dig,st(i,1),100,0,ceil(trial_end(i)-st(i,1))*1000); % sampling
    raw{i,3}=c.spikeraster.times;
    c=spike2eventRasteandPSTH_NP(data.Up.dig_offset,st(i,1),100,0,ceil(trial_end(i)-st(i,1))*1000); % timestamps of Decision onset
    raw{i,4}=c.spikeraster.times;
    c=spike2eventRasteandPSTH_NP(data.LeftSp.dig,st(i,1),100,0,ceil(trial_end(i)-st(i,1))*1000); % timestamps of Left licks
    raw{i,5}=c.spikeraster.times;
    c=spike2eventRasteandPSTH_NP(data.RightSp.dig,st(i,1),100,0,ceil(trial_end(i)-st(i,1))*1000); % timestamps of right licks
    raw{i,6}=c.spikeraster.times;
    if isempty(raw{i,5}) & isempty(raw{i,6})
        trial_noResponse=[trial_noResponse i];
    end
    switch raw{i,1}
        case 1 % left choice
            if isempty(raw{i,5}) & ~isempty(raw{i,6})
                LeftError.trial=[LeftError.trial i];
            end
            if ~isempty(raw{i,5}) & ~isempty(raw{i,6})
                if raw{i,6}(1) < raw{i,5}(1)
                    LeftError.trial=[LeftError.trial i];
                end
            end
        case 2 % right choice
            if isempty(raw{i,6}) & ~isempty(raw{i,5})
                RightError.trial=[RightError.trial i];
            end
            if ~isempty(raw{i,5}) & ~isempty(raw{i,6})
                if raw{i,5}(1) < raw{i,6}(1)
                    RightError.trial=[RightError.trial i];
                end
            end
    end
    raw{i,7}=st(i,3);
end

LeftError.num=length(LeftError.trial);
RightError.num=length(RightError.trial);

%% 2nd find out trial with no decision making
sumData.noResponse=trial_noResponse;
left_trial=find(st(:,2)==1);
left_trial_decision=setdiff(left_trial,trial_noResponse);
Right_trial=find(st(:,2)==2);
Right_trial_decision=setdiff(Right_trial,trial_noResponse);
%% 3.2 Let's extract the "correct response" with first lick
trial_decision=[left_trial_decision; Right_trial_decision];
trial_decision=sort(trial_decision);
trial_decision(:,2)=st(trial_decision,2);
trial_decision(:,3)=1;
% sum all error trial
if isfield(LeftError,'trial') & isfield(RightError,'trial')
    error_all=[LeftError.trial RightError.trial];
elseif isfield(LeftError,'trial')
    error_all=LeftError.trial;
elseif isfield(RightError,'trial')
    error_all=RightError.trial;
else
    error_all=[];
end
if isempty(error_all)
else
    error_all=sort(error_all);
    for i=1:length(error_all)
        index=find(trial_decision(:,1)==error_all(i));
        trial_decision(index,3)=0;
    end
end
B=cumsum(trial_decision(:,3));
for i=1:size(trial_decision,1)
    trial_decision(i,4)=B(i)/i;
end
%% save data
sumData.LeftTrialDecision=left_trial_decision;
sumData.LeftError=LeftError;
sumData.RightTrialDecision=Right_trial_decision;
sumData.RightError=RightError;
sumData.MissLateral=length(trial_noResponse)/(size(st,1));
sumData.Performance=1-(LeftError.num+RightError.num)/length([left_trial_decision;Right_trial_decision]);
sumData.Performance_bias=LeftError.num/length(left_trial_decision)-RightError.num/length(Right_trial_decision);
sumData.Trial_decision=trial_decision;
sumData.raw=raw;
sumData.raw_decision=raw(trial_decision(:,1),:);
% clearvars -except sumData FinalPerformance data tasteID x cwd
% save('sumData.mat','sumData')
a=num2cell(sumData.Trial_decision);
sumData.raw_decision=[a sumData.raw_decision];
sumData.raw_decision(:,5)=[];
sumData.tasteID          = x;
%% load optogenetic stimulation information
thr = 2;
file = dir('*.rhd');
dataRaw = read_Intan(file.name);
[data.OptSt_on,data.OptSt_off] = Timing_onset_offset(dataRaw.analog(4,:), dataRaw.ts, thr,300,0);
trial = sumData;
%%
clear Left Right
for i = 1: size(trial.raw_decision,1)
    if ~isempty(trial.raw_decision{i,6})
        central_last(i) = trial.raw_decision{i,6}(end) + trial.raw_decision{i,5};
        central_first(i) = trial.raw_decision{i,6}(1) + trial.raw_decision{i,5};
        Left =  trial.raw_decision{i,8};
        Right = trial.raw_decision{i,9};
        lateral(i) = min([Left,Right])+ trial.raw_decision{i,5};
    end
end
% sampling_duration = central_last - central_first;
% delay_duration    = lateral-central_last;
k=1;
%%
switch delay_idx % whether it is for delay epoch.
    %% for delay
    case 1
        for i = 1: size(trial.raw_decision,1)
            for j = 1:length(data.OptSt_on)
                if data.OptSt_on(j)> central_last(i)-0.2 && data.OptSt_on(j)< central_last(i)+2 % edit by Ke to add -0.2 to the central last
                    sti(k) = i;
                    k = k+1;
                end
            end
        end
    case 0
        %% for detection
        for i = 1: size(trial.raw_decision,1)
            for j = 1:length(data.OptSt_on)
                if data.OptSt_on(j)< central_first(i) && data.OptSt_on(j)> central_first(i)-2
                    sti(k) = i;
                    k = k+1;
                end
            end
        end
end
sumData.stiTrial = sti;
performance_sti(:,1) = cell2mat(trial.raw_decision(sti,3));
performance_sti(:,2) = cell2mat(trial.raw_decision(sti,2));
L_R = length(find(performance_sti(:,2)==1))/(length(find(performance_sti(:,2)==1))+length(find(performance_sti(:,2)==2)));
avg_per_sti = sum(performance_sti(:,1))/size(performance_sti,1);
%% get the performance of the trials without stimulation
perf = cell2mat(trial.raw_decision(:,3));
perf(sti)=[];
avg_per_non = sum(perf)/length(perf);
sumData.avg_per_sti = avg_per_sti; % trials with light stimulation
sumData.avg_per_non =avg_per_non;  % trials without light stimulation
