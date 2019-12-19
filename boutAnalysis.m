%% Mannually import the interlick interval data
% Store the interlick interval from davis into the data
data = NH;
%%
for i = 1:size(data,1)
    info(i) = bout(data(i,:));
end
%%
for i = 1:size(data,1)
    info(i).boutnum = length(info(i).duration);
end
clearvars -except info
%% Manually import trial information to sort trials
% Store the trial information into trial
trial = NH(:,3);
trialtype = unique(trial);

for i = 1:length(trialtype)
    trialSe{1,i} = find(trial == trialtype(i));
end
infozero =[];
for i = 1:length(trialSe)
    idx = trialSe{1,i};
    for j = 1:length(idx)
        infozero = [infozero,info(idx(j))];
    end
    trialSe{2,i}= infozero;
    infozero =[];
end

%%