function data = iti_decision(data)
for i = 1:length(data)
     temp = data(i).raw_decision; % extract the raw decision information
     for j = 1:size(temp,1)
         iti(j) = mean(diff(temp{j,6})); % calculate the sampling frequency for each trial
     end
%      iti(isnan(iti))=[];
     tasteID = unique(cell2mat(temp(:,10))); % find the tasteID line
     taste = data(i).tasteID; % find the tasteID
     for z = 1:length(taste)
         temp2 = iti(find(cell2mat(temp(:,10)) == tasteID(z)));  
         temp2(isnan(temp2)) = [];
         data(i).([taste{z},'_iti']) = mean(temp2);
     end
%      data(i).sampling = mean(sampling);
     clear iti temp2
end