function data = iti_decision_lateral(data)
tr = {'left', 'right'}
for i = 1:length(data)
     temp = data(i).raw_decision; % extract the raw decision information
     temp2 = cell2mat(temp(:,2:3));
     for j = 1:length(tr)
         idx.(tr{j}) = find(temp2(:,1) ==j & temp2(:,2) ==1);
         for z = 1:length(idx.(tr{j}))
             lateral = temp{idx.(tr{j})(z),j+7};
             central = temp{idx.(tr{j})(z),6};
%              delay.(tr{j})(z) = lateral(1) - central(end);
             decision.(tr{j})(z) = mean(diff(lateral));
         end
%          data(i).([tr{j},'delay']) = mean([delay.(tr{j})]);
         data(i).([tr{j},'decision_iti']) = mean([decision.(tr{j})],'omitnan');
     end
     clear decision
end
clearvars -except data 