%%
load('DataForLicking.mat')
%% calculate the sampling duration for each tastants
clearvars -except data
for i = 1:length(data)
     temp = data(i).raw_decision; % extract the raw decision information
     for j = 1:size(temp,1)
         sampling(j) = temp{j,6}(end) - temp{j,6}(1); % calculate the sampling duration for each trial
     end
     tasteID = unique(cell2mat(temp(:,10))); % find the tasteID line
     taste = data(i).tasteID; % find the tasteID
     for z = 1:length(taste)
         data(i).(taste{z}) = mean(sampling(find(cell2mat(temp(:,10)) == tasteID(z))));  
     end
     data(i).sampling = mean(sampling);
     clear sampling
end

data = iti_decision(data);

%%%% summary the sampling duration and InterLick Interval for each tastant
for i = 1:length(data)
    id{i} = data(i).mouseID;
end
id_uni = unique(id);
% sd_f = {'S', 'M', 'SO', 'Q','sampling'}
sd_f = {'S', 'M', 'SO', 'Q'};

for i = 1:length(id_uni)
    duration(i).id = id_uni{i};
    iti(i).id      = id_uni{i};
    idx = find(strcmp(id,id_uni{i}));
    if length(idx) >=1
        for j = 1:length(sd_f)
            duration(i).(sd_f{j}) = mean([data(idx).(sd_f{j})]);
            iti(i).(sd_f{j})      = mean([data(idx).([sd_f{j},'_iti'])]);
        end
    else
        error('Something is wrong')
    end
end
clearvars -except data duration

%% Save the results as table 
iti = struct2table(iti);
duration_temp = struct2table(duration);

writetable(duration_temp,'Central-licking.xls','sheet',1)
writetable(iti,'Central-licking.xls','sheet',2)

%% plot the result and perform the statistic test
temp = duration;
temp = rmfield(temp,'id');
temp = cell2mat(squeeze(struct2cell(temp))');

control = temp;
bar_plot_multi(control)
name = {[],'S','M','SO','Q',[]}
set(gca,'XTick',0:5)
set(gca,'xticklabel',name)
ylabel('Sampling duration (s)')

%%%%% perform one-way anova
% temp = temp(:,1:4); % use compare duration for each tastant
% p = anova1(temp)
%% Calculate the delay 
close all
clearvars -except data duration
% Only use the correct trials to calculate the delay and lateral licking
% duration
tr = {'left', 'right'}
for i = 1:length(data)
     temp = data(i).raw_decision; % extract the raw decision information
     temp2 = cell2mat(temp(:,2:3));
     for j = 1:length(tr)
         idx.(tr{j}) = find(temp2(:,1) ==j & temp2(:,2) ==1);
         for z = 1:length(idx.(tr{j}))
             lateral = temp{idx.(tr{j})(z),j+7};
             central = temp{idx.(tr{j})(z),6};
             delay.(tr{j})(z) = lateral(1) - central(end);
             decision.(tr{j})(z) = lateral(end)-lateral(1);
         end
         data(i).([tr{j},'delay']) = mean([delay.(tr{j})]);
         data(i).([tr{j},'decision']) = mean([decision.(tr{j})]);
     end
     clear delay decision
end
clearvars -except data duration
% summary the delay and decision for left/right trials
for i = 1:length(data)
    id{i} = data(i).mouseID;
end
id_uni = unique(id);
sd_f = {'leftdelay', 'leftdecision', 'rightdelay', 'rightdecision'};
for i = 1:length(id_uni)
    timing(i).id = id_uni{i};
    idx = find(strcmp(id,id_uni{i}));
    if length(idx) >=1
        for j = 1:length(sd_f)
            timing(i).(sd_f{j}) = mean([data(idx).(sd_f{j})]);
        end
    else
        error('Something is wrong')
    end
end
clearvars -except data duration timing

%% plot the duration for the delay
temp = timing;
temp = rmfield(temp,'id');
temp = cell2mat(squeeze(struct2cell(temp))'); % the first collumn is the left delay and the 3rd collumn is the right delay

control = temp;
barplot_test(control(:,1), control(:,3))
name = {[],'Left','Right',[]}
set(gca,'xticklabel',name)
ylabel('Delay (s)')
ylim([0,3])
% perform the t-test
% [h,p] = ttest2(temp(:,1),temp(:,3));
%% plot the duration for lateral licking
barplot_test(control(:,2), control(:,4))
set(gca,'xticklabel',name)
ylabel('Lateral licking (s)')
ylim([0,2])
% perform the t-test
% [h,p] = ttest2(temp(:,2),temp(:,4));

%% Save the results as table 
timing = struct2table(timing);
timing_temp = struct2table(duration);

writetable(timing,'Timing_delay_lateral.xls','sheet',1)
% writetable(iti,'Central-licking.xls','sheet',2)
%% calculate the performance for each tastant
clearvars -except data duration timing
for i = 1:length(data)
     temp = data(i).raw_decision; % extract the raw decision information
     performance = cell2mat(temp(:,3));
     tasteID = unique(cell2mat(temp(:,10))); % find the tasteID line
     taste = data(i).tasteID; % find the tasteID
     for z = 1:length(taste)
         t = cell2mat(temp(:,10));
         p = performance(find(t == tasteID(z)));
         data(i).([(taste{z}),'_p']) = sum(p)/length(p);  
     end
end


%%%% summary the sampling duration for each tastant
clear performance
for i = 1:length(data)
    id{i} = data(i).mouseID;
end
id_uni = unique(id);
sd_f = {'S_p', 'M_p', 'SO_p', 'Q_p'}
for i = 1:length(id_uni)
    performance(i).id = id_uni{i};
    idx = find(strcmp(id,id_uni{i}));
    if length(idx) >=1
        for j = 1:length(sd_f)
            performance(i).(sd_f{j}) = mean([data(idx).(sd_f{j})]);
        end
    else
        error('Something is wrong')
    end
end
clearvars -except data duration timing performance


%% plot the result and perform the statistic test
temp = performance;
temp = rmfield(temp,'id');
temp = cell2mat(squeeze(struct2cell(temp))');
bar_plot_multi(temp)
name = {[],'S','M','SO','Q',[]}
set(gca,'xticklabel',name)
ylabel('Performance')
ylim([0,1])
%% perform one-way anova
temp = temp(:,1:4); % use compare duration for each tastant
p = anova1(temp)
%%timing = struct2table(timing);
performance = struct2table(performance);
writetable(performance,'Performance.xls','sheet',1)