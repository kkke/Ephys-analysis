%%
load('DataForLicking.mat')
%% calculate the sampling duration for each tastants
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


%%%% summary the sampling duration for each tastant
for i = 1:length(data)
    id{i} = data(i).mouseID;
end
id_uni = unique(id);
sd_f = {'S', 'M', 'SO', 'Q','sampling'}
for i = 1:length(id_uni)
    duration(i).id = id_uni{i};
    idx = find(strcmp(id,id_uni{i}));
    if length(idx) >=1
        for j = 1:length(sd_f)
            duration(i).(sd_f{j}) = mean([data(idx).(sd_f{j})]);
        end
    else
        error('Something is wrong')
    end
end
clearvars -except data duration


%%% plot the result and perform the statistic test
temp = duration;
temp = rmfield(temp,'id');
temp = cell2mat(squeeze(struct2cell(temp))');

control = temp;
m1 = mean(control(:,1));
m2 = mean(control(:,2));
m3 = mean(control(:,3));
m4 = mean(control(:,4));

name = {[],'S','M','SO','Q',[]}
sem1 = std(control(:,1))./sqrt(length(control(:,3)));
sem2 = std(control(:,2))./sqrt(length(control(:,3)));
sem3 = std(control(:,3))./sqrt(length(control(:,3)));
sem4 = std(control(:,4))./sqrt(length(control(:,4)));

figure;
bar(1,m1,'EdgeColor',[0 0 0],'FaceColor',[0.5,0.5,1],'LineWidth',1);hold on
bar(2,m2,'EdgeColor',[0 0 0],'FaceColor',[1,0.5,0.25],'LineWidth',1);hold on
bar(3,m3,'EdgeColor',[0 0 0],'FaceColor',[1,0.5,0.25],'LineWidth',1);hold on
bar(4,m4,'EdgeColor',[0 0 0],'FaceColor',[1,0.5,0.25],'LineWidth',1);hold on

set(gca,'XTick',0:5)
set(gca,'xticklabel',name)
ylabel('Sampling duration (s)')
hold on;
e=errorbar([m1, m2, m3, m4], [sem1, sem2, sem3, sem4],'LineWidth',1.5);
% plot([1 1],[mean(a_m) mean(a_m)+a_sem],'r');hold on;
e.Color ='k';
e.LineStyle = 'none';

%%%%% perform one-way anova
temp = temp(:,1:4); % use compare duration for each tastant
p = anova1(temp)
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


% plot the result and perform the statistic test
temp = timing;
temp = rmfield(temp,'id');
temp = cell2mat(squeeze(struct2cell(temp))'); % the first collumn is the left delay and the 3rd collumn is the right delay

control = temp;
m1 = mean(control(:,1));
m2 = mean(control(:,3));
name = {[],'Left','Right',[]}
sem1 = std(control(:,1))./sqrt(length(control(:,1)));
sem2 = std(control(:,3))./sqrt(length(control(:,3)));
figure;
bar(1,m1,'EdgeColor',[0 0 0],'FaceColor',[0.5,0.5,1],'LineWidth',1);hold on
bar(2,m2,'EdgeColor',[0 0 0],'FaceColor',[1,0.5,0.25],'LineWidth',1);hold on
set(gca,'XTick',0:3)
set(gca,'xticklabel',name)
ylabel('Delay (s)')
hold on;
e=errorbar([m1, m2], [sem1, sem2],'LineWidth',1.5);
e.Color ='k';
e.LineStyle = 'none';
ylim([0,3])
% perform the t-test
[h,p] = ttest2(temp(:,1),temp(:,3));
%%
figure;
m3 = mean(control(:,2));
m4 = mean(control(:,4));
sem3 = std(control(:,2))./sqrt(length(control(:,2)));
sem4 = std(control(:,4))./sqrt(length(control(:,4)));
bar(1,m3,'EdgeColor',[0 0 0],'FaceColor',[1,0.5,0.25],'LineWidth',1);hold on
bar(2,m4,'EdgeColor',[0 0 0],'FaceColor',[1,0.5,0.25],'LineWidth',1);hold on
set(gca,'XTick',0:3)
set(gca,'xticklabel',name)
ylabel('Lateral licking (s)')
hold on;
e=errorbar([m3, m4], [ sem3, sem4],'LineWidth',1.5);
% plot([1 1],[mean(a_m) mean(a_m)+a_sem],'r');hold on;
e.Color ='k';
e.LineStyle = 'none';

% perform the t-test
[h,p] = ttest2(temp(:,2),temp(:,4));
ylim([0,2])
%%
clearvars -except data duration timing
save('DataForLicking.mat','data','duration','timing')