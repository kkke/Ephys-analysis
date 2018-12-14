%% Trying extract the performance for each individual tastant
% Let take RVKC205 as an example
% S is line 1, M is line 3, Q is line 4, SO is line 5
clear p
taste = {'M','S','Q','SO'}
id    = [2,3,4,5];
for i = 1:length(summaryData)
    tasteId = cell2mat(summaryData(i).raw_decision(:,10));
    for j = 1:length(taste)
        per = cell2mat(summaryData(i).raw_decision(find(tasteId == id(j)),3));
        p(i).(taste{j}) = sum(per)/length(per);
    end   
end
%% save performance for each taste
save('Performance.mat','p')
