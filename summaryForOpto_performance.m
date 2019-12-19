taste = {'S', 'Q', 'M','SO'};

%% get the stimulation trial
for i = 1:length(data_delay)
    for j = 1:length(taste)
        temp = data_delay(i).raw_decision(data_delay(i).stiTrial,:);
%         temp= data_delay(i).raw_decision;
%         temp(data_delay(i).stiTrial,:)=[];
        idx = find(strcmp(data_delay(1).tasteID,taste{j}));
        temp_idx = find(cell2mat(temp(:,end))==idx);
        performance(i).(['sti_',(taste{j})])= sum(cell2mat(temp(temp_idx,3)))/length(cell2mat(temp(temp_idx,3)));      
    end
    clear temp idx temp_idx
end

for i = 1:length(data_delay)
    performance(i).sti_sweet  = (performance(i).sti_S + performance(i).sti_M)/2;
    performance(i).sti_bitter = (performance(i).sti_Q + performance(i).sti_SO)/2;
end
%% get the non-stimulation trial
for i = 1:length(data_delay)
    for j = 1:length(taste)
%         temp = data_delay(i).raw_decision(data_delay(i).stiTrial,:);
        temp= data_delay(i).raw_decision;
        temp(data_delay(i).stiTrial,:)=[];
        idx = find(strcmp(data_delay(1).tasteID,taste{j}));
        temp_idx = find(cell2mat(temp(:,end))==idx);
        performance(i).(['nonsti_',(taste{j})])= sum(cell2mat(temp(temp_idx,3)))/length(cell2mat(temp(temp_idx,3)));      
    end
    clear temp idx temp_idx
end

for i = 1:length(data_delay)
    performance(i).nonsti_sweet  = (performance(i).nonsti_S + performance(i).nonsti_M)/2;
    performance(i).nonsti_bitter = (performance(i).nonsti_Q + performance(i).nonsti_SO)/2;
end

for i = 1:length(data_delay)
    performance(i).relative_sweet = 1-performance(i).sti_sweet/performance(i).nonsti_sweet;
    performance(i).relative_bitter = 1-performance(i).sti_bitter/performance(i).nonsti_bitter;
end