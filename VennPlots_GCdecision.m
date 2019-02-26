H=1;
for j=1:214
    if tasteFlag4(1,j)==1 & sum(tasteFlag2(j,:))>0
        summary(j).taste=1;
    else
        summary(j).taste=0;
    end
end
%%
clearvars -except summary
idx = intersect(find([summary.planning]==1), find([summary.action]==1));

idxx = intersect(idx, find([summary.taste]==1));

idxxx = intersect(find([summary.planning]==1), find([summary.taste]==1));

idxxxx = intersect(find([summary.action]==1), find([summary.taste]==1));

sum([summary.planning])
sum([summary.action])
sum([summary.taste])
sum([summary.outcome])

%%
clearvars -except summary
for i = 1:length(summary)
    summary(i).tastePlan = 0;
end
plan = find([summary.planning]==1);
for i = 1:length(bothnon)
    
    summary(plan(bothnon(i))).tastePlan = 1;
end
%%
clearvars -except summary
idx = intersect(find([summary.planning]==1), find([summary.tastePlan]==1));

idxx = intersect(idx, find([summary.taste]==1));

%%
clearvars -except summary
idx = intersect(find([summary.action]==1), find([summary.taste]==1));
idxx = intersect(idx, find([summary.outcome]==1));
idxxx = intersect(find([summary.action]==1), find([summary.outcome]==1));
idxxxx = intersect(find([summary.taste]==1), find([summary.outcome]==1));
