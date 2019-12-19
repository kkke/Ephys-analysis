%% Trying to calculate the outcome response
% first let's get the event for correct and error trials
for i = 1:length(Sum)
    outcome(i).corr    = [Sum(i).event.tsRCorr.FLickRSpou; Sum(i).event.tsLCorr.FLickLSpou];
    outcome(i).error   = [Sum(i).event.tsLErr.FLickRSpou; Sum(i).event.tsRErr.FLickRSpou];
    [outcome(i).dscorr, outcome(i).idx]  = datasample(outcome(i).corr,length(outcome(i).error),'Replace',false);
    temp = [];
    for j = 1:size(Sum(i).event.raw_decision,1)
        if Sum(i).event.raw_decision{j,3} == 1 & Sum(i).event.raw_decision{j,2} ==1
            temp{j,1} =  Sum(i).event.raw_decision{j,5} + Sum(i).event.raw_decision{j,8}(1);
            temp{j,5} =  (Sum(i).event.raw_decision{j,5} + Sum(i).event.raw_decision{j,8})';
        elseif Sum(i).event.raw_decision{j,3} == 1 & Sum(i).event.raw_decision{j,2} ==2
            temp{j,1} =  Sum(i).event.raw_decision{j,5} + Sum(i).event.raw_decision{j,9}(1);
            temp{j,5} =  (Sum(i).event.raw_decision{j,5} + Sum(i).event.raw_decision{j,9})';
        elseif Sum(i).event.raw_decision{j,3} == 0 & Sum(i).event.raw_decision{j,2} ==1
            temp{j,1} =  Sum(i).event.raw_decision{j,5} + Sum(i).event.raw_decision{j,9}(1);
            temp{j,6} =  (Sum(i).event.raw_decision{j,5} + Sum(i).event.raw_decision{j,9})';

        else   Sum(i).event.raw_decision{j,3} == 0 & Sum(i).event.raw_decision{j,2} ==2
            temp{j,1} =  Sum(i).event.raw_decision{j,5} + Sum(i).event.raw_decision{j,8}(1);
            temp{j,6} =  (Sum(i).event.raw_decision{j,5} + Sum(i).event.raw_decision{j,8})';

        end
        temp{j,2} =  Sum(i).event.raw_decision{j,5} + Sum(i).event.raw_decision{j,8};
        temp{j,3} =  Sum(i).event.raw_decision{j,5} + Sum(i).event.raw_decision{j,9};
        temp{j,4} =  Sum(i).event.raw_decision{j,3};
    end
    outcome(i).temp = temp;
end
%% find the licks for correct and error trials
for i = 1:length(outcome)
    corLicks = outcome(i).temp(:,5);
    [~,~,ib] = intersect(outcome(i).dscorr, cell2mat(outcome(i).temp(:,1)));
    outcome(i).corLicks = cell2mat(corLicks(ib));
    errorib  = find(cell2mat(outcome(i).temp(:,4)) == 0);
    errLicks = outcome(i).temp(:,6);
    outcome(i).errLicks = cell2mat(errLicks(errorib));
end
%% lets compute the perference
for i = 1:length(Sum)
    fprintf('Process Neurons %4.2f \n', i)
    unit.corr   = spike2eventRasteandPSTH_NP(Sum(i).timestampN,outcome(i).dscorr,500,0,500);
    unit.error  = spike2eventRasteandPSTH_NP(Sum(i).timestampN,outcome(i).error,500,0,500);
    p(i) = Test_auROC_dISCRIMINATION(unit.corr.scmatrix,unit.error.scmatrix);
    % Using bootstrap to perform the statistic test
    com_data = [unit.corr .scmatrix; unit.error.scmatrix];
    for j = 1:1000 % 1000 iteration
        [data(:,j),idx_p] = datasample(com_data,length(unit.corr.scmatrix),'Replace',false);
        idx_pall          = 1:length(com_data);
        data_left(:,j)    = com_data(setdiff(idx_pall,idx_p));
        outcome(i).stats.shaffle_p(j)      = Test_auROC_dISCRIMINATION(data(:,j),data_left(:,j));
    end
    if p(i)>=0 % one tail bootstrap
        outcome(i).stats.p = length(find(outcome(i).stats.shaffle_p>=p(i)))/1000;
    else
        outcome(i).stats.p = length(find(outcome(i).stats.shaffle_p<=p(i)))/1000;
    end
    
    if outcome(i).stats.p>0.01
        outcome(i).stats.resp = 0;
    else
        outcome(i).stats.resp = 1;
    end
    clear data data_left shaffle_p
end
%% let's see how many neurons encode outcome
for i = 1:length(outcome)
    outcome(i).resp = outcome(i).stats.resp;
end
%% 
%% lets compute the perference just with licking data
for i = 1:length(Sum)
    fprintf('Process Neurons %4.2f \n', i)
    unit.corr   = spike2eventRasteandPSTH_NP(outcome(i).corLicks,outcome(i).dscorr,500,0,500);
    unit.error  = spike2eventRasteandPSTH_NP(outcome(i).errLicks,outcome(i).error,500,0,500);
    p(i) = Test_auROC_dISCRIMINATION(unit.corr.scmatrix,unit.error.scmatrix);
    % Using bootstrap to perform the statistic test
    com_data = [unit.corr .scmatrix; unit.error.scmatrix];
    for j = 1:1000 % 1000 iteration
        [data(:,j),idx_p] = datasample(com_data,length(unit.corr.scmatrix),'Replace',false);
        idx_pall          = 1:length(com_data);
        data_left(:,j)    = com_data(setdiff(idx_pall,idx_p));
        outcome(i).Lickstats.shaffle_p(j)      = Test_auROC_dISCRIMINATION(data(:,j),data_left(:,j));
    end
    if p(i)>=0 % one tail bootstrap
        outcome(i).Lickstats.p = length(find(outcome(i).Lickstats.shaffle_p>=p(i)))/1000;
    else
        outcome(i).Lickstats.p = length(find(outcome(i).Lickstats.shaffle_p<=p(i)))/1000;
    end
    
    if outcome(i).Lickstats.p>0.01
        outcome(i).Lickstats.resp = 0;
    else
        outcome(i).Lickstats.resp = 1;
    end
    clear data data_left shaffle_p
end
%% try to use the pesudo correct and error licks
temp =[];
for i = 1:length(outcome)
    for j = 1:size(outcome(i).temp,1)
        temp{j,1} = ([outcome(i).temp{j,2},outcome(i).temp{j,3}])';
    end
%     Licks = outcome(i).temp(:,5);
    [~,~,ib] = intersect(outcome(i).dscorr, cell2mat(outcome(i).temp(:,1)));
    outcome(i).pseucorLicks = cell2mat(temp(ib));
    errorib  = find(cell2mat(outcome(i).temp(:,4)) == 0);
%     errLicks = outcome(i).temp(:,6);
    outcome(i).pseuerrLicks = cell2mat(temp(errorib));
end
%% lets compute the perference just with licking data
% p means 500 ms
% pp means 200 ms
% ppp means 100 ms
% pppp means 300 ms
% ppppp means 400 ms
clear data data_left shaffle_p p

for i = 1:length(Sum)
    fprintf('Process Neurons %4.2f \n', i)
    unit.corr   = spike2eventRasteandPSTH_NP(outcome(i).pseucorLicks,outcome(i).dscorr,200,0,200);
    unit.error  = spike2eventRasteandPSTH_NP(outcome(i).pseuerrLicks,outcome(i).error,200,0,200);
    p(i) = Test_auROC_dISCRIMINATION(unit.corr.scmatrix,unit.error.scmatrix);
    % Using bootstrap to perform the statistic test
    com_data = [unit.corr .scmatrix; unit.error.scmatrix];
    for j = 1:1000 % 1000 iteration
        [data(:,j),idx_p] = datasample(com_data,length(unit.corr.scmatrix),'Replace',false);
        idx_pall          = 1:length(com_data);
        data_left(:,j)    = com_data(setdiff(idx_pall,idx_p));
        outcome(i).ppLickstats.shaffle_p(j)      = Test_auROC_dISCRIMINATION(data(:,j),data_left(:,j));
    end
    if p(i)>=0 % one tail bootstrap
        outcome(i).ppLickstats.p = length(find(outcome(i).ppLickstats.shaffle_p>=p(i)))/1000;
    else
        outcome(i).ppLickstats.p = length(find(outcome(i).ppLickstats.shaffle_p<=p(i)))/1000;
    end
    
    if outcome(i).ppLickstats.p>0.01
        outcome(i).ppLickstats.resp = 0;
    else
        outcome(i).ppLickstats.resp = 1;
    end
    clear data data_left shaffle_p 
end
%
for i = 1:length(Sum)
    fprintf('Process Neurons %4.2f \n', i)
    unit.corr   = spike2eventRasteandPSTH_NP(Sum(i).timestampN,outcome(i).dscorr,200,0,200);
    unit.error  = spike2eventRasteandPSTH_NP(Sum(i).timestampN,outcome(i).error,200,0,200);
    p(i) = Test_auROC_dISCRIMINATION(unit.corr.scmatrix,unit.error.scmatrix);
    % Using bootstrap to perform the statistic test
    com_data = [unit.corr .scmatrix; unit.error.scmatrix];
    for j = 1:1000 % 1000 iteration
        [data(:,j),idx_p] = datasample(com_data,length(unit.corr.scmatrix),'Replace',false);
        idx_pall          = 1:length(com_data);
        data_left(:,j)    = com_data(setdiff(idx_pall,idx_p));
        outcome(i).ppLickstats.shaffle_p(j)      = Test_auROC_dISCRIMINATION(data(:,j),data_left(:,j));
    end
    if p(i)>=0 % one tail bootstrap
        outcome(i).ppstats.p = length(find(outcome(i).ppLickstats.shaffle_p>=p(i)))/1000;
    else
        outcome(i).ppstats.p = length(find(outcome(i).ppLickstats.shaffle_p<=p(i)))/1000;
    end
    
    if outcome(i).ppstats.p>0.01
        outcome(i).ppstats.resp = 0;
    else
        outcome(i).ppstats.resp = 1;
    end
    clear data data_left shaffle_p 
end
for i = 1:length(outcome)
    outcome(i).ppresp = outcome(i).ppstats.resp;
end
%%
clear data data_left shaffle_p p

for i = 1:length(Sum)
    fprintf('Process Neurons %4.2f \n', i)
    unit.corr   = spike2eventRasteandPSTH_NP(outcome(i).pseucorLicks,outcome(i).dscorr,100,0,100);
    unit.error  = spike2eventRasteandPSTH_NP(outcome(i).pseuerrLicks,outcome(i).error,100,0,100);
    p(i) = Test_auROC_dISCRIMINATION(unit.corr.scmatrix,unit.error.scmatrix);
    % Using bootstrap to perform the statistic test
    com_data = [unit.corr .scmatrix; unit.error.scmatrix];
    for j = 1:1000 % 1000 iteration
        [data(:,j),idx_p] = datasample(com_data,length(unit.corr.scmatrix),'Replace',false);
        idx_pall          = 1:length(com_data);
        data_left(:,j)    = com_data(setdiff(idx_pall,idx_p));
        outcome(i).pppLickstats.shaffle_p(j)      = Test_auROC_dISCRIMINATION(data(:,j),data_left(:,j));
    end
    if p(i)>=0 % one tail bootstrap
        outcome(i).pppLickstats.p = length(find(outcome(i).pppLickstats.shaffle_p>=p(i)))/1000;
    else
        outcome(i).pppLickstats.p = length(find(outcome(i).pppLickstats.shaffle_p<=p(i)))/1000;
    end
    
    if outcome(i).pppLickstats.p>0.01
        outcome(i).pppLickstats.resp = 0;
    else
        outcome(i).pppLickstats.resp = 1;
    end
    clear data data_left shaffle_p 
end
% clear p
for i = 1:length(Sum)
    fprintf('Process Neurons %4.2f \n', i)
    unit.corr   = spike2eventRasteandPSTH_NP(Sum(i).timestampN,outcome(i).dscorr,100,0,100);
    unit.error  = spike2eventRasteandPSTH_NP(Sum(i).timestampN,outcome(i).error,100,0,100);
    p(i) = Test_auROC_dISCRIMINATION(unit.corr.scmatrix,unit.error.scmatrix);
    % Using bootstrap to perform the statistic test
    com_data = [unit.corr .scmatrix; unit.error.scmatrix];
    for j = 1:1000 % 1000 iteration
        [data(:,j),idx_p] = datasample(com_data,length(unit.corr.scmatrix),'Replace',false);
        idx_pall          = 1:length(com_data);
        data_left(:,j)    = com_data(setdiff(idx_pall,idx_p));
        outcome(i).pppLickstats.shaffle_p(j)      = Test_auROC_dISCRIMINATION(data(:,j),data_left(:,j));
    end
    if p(i)>=0 % one tail bootstrap
        outcome(i).pppstats.p = length(find(outcome(i).pppLickstats.shaffle_p>=p(i)))/1000;
    else
        outcome(i).pppstats.p = length(find(outcome(i).pppLickstats.shaffle_p<=p(i)))/1000;
    end
    
    if outcome(i).pppstats.p>0.01
        outcome(i).pppstats.resp = 0;
    else
        outcome(i).pppstats.resp = 1;
    end
    clear data data_left shaffle_p 
end
for i = 1:length(outcome)
    outcome(i).pppresp = outcome(i).pppstats.resp;
end
%%
clear data data_left shaffle_p p

for i = 1:length(Sum)
    fprintf('Process Neurons %4.2f \n', i)
    unit.corr   = spike2eventRasteandPSTH_NP(outcome(i).pseucorLicks,outcome(i).dscorr,300,0,300);
    unit.error  = spike2eventRasteandPSTH_NP(outcome(i).pseuerrLicks,outcome(i).error,300,0,300);
    p(i) = Test_auROC_dISCRIMINATION(unit.corr.scmatrix,unit.error.scmatrix);
    % Using bootstrap to perform the statistic test
    com_data = [unit.corr .scmatrix; unit.error.scmatrix];
    for j = 1:1000 % 1000 iteration
        [data(:,j),idx_p] = datasample(com_data,length(unit.corr.scmatrix),'Replace',false);
        idx_pall          = 1:length(com_data);
        data_left(:,j)    = com_data(setdiff(idx_pall,idx_p));
        outcome(i).ppppLickstats.shaffle_p(j)      = Test_auROC_dISCRIMINATION(data(:,j),data_left(:,j));
    end
    if p(i)>=0 % one tail bootstrap
        outcome(i).ppppLickstats.p = length(find(outcome(i).ppppLickstats.shaffle_p>=p(i)))/1000;
    else
        outcome(i).ppppLickstats.p = length(find(outcome(i).ppppLickstats.shaffle_p<=p(i)))/1000;
    end
    
    if outcome(i).ppppLickstats.p>0.01
        outcome(i).ppppLickstats.resp = 0;
    else
        outcome(i).ppppLickstats.resp = 1;
    end
    clear data data_left shaffle_p 
end
% clear p
for i = 1:length(Sum)
    fprintf('Process Neurons %4.2f \n', i)
    unit.corr   = spike2eventRasteandPSTH_NP(Sum(i).timestampN,outcome(i).dscorr,300,0,300);
    unit.error  = spike2eventRasteandPSTH_NP(Sum(i).timestampN,outcome(i).error,300,0,300);
    p(i) = Test_auROC_dISCRIMINATION(unit.corr.scmatrix,unit.error.scmatrix);
    % Using bootstrap to perform the statistic test
    com_data = [unit.corr .scmatrix; unit.error.scmatrix];
    for j = 1:1000 % 1000 iteration
        [data(:,j),idx_p] = datasample(com_data,length(unit.corr.scmatrix),'Replace',false);
        idx_pall          = 1:length(com_data);
        data_left(:,j)    = com_data(setdiff(idx_pall,idx_p));
        outcome(i).ppppLickstats.shaffle_p(j)      = Test_auROC_dISCRIMINATION(data(:,j),data_left(:,j));
    end
    if p(i)>=0 % one tail bootstrap
        outcome(i).ppppstats.p = length(find(outcome(i).ppppLickstats.shaffle_p>=p(i)))/1000;
    else
        outcome(i).ppppstats.p = length(find(outcome(i).ppppLickstats.shaffle_p<=p(i)))/1000;
    end
    
    if outcome(i).ppppstats.p>0.01
        outcome(i).ppppstats.resp = 0;
    else
        outcome(i).ppppstats.resp = 1;
    end
    clear data data_left shaffle_p 
end
for i = 1:length(outcome)
    outcome(i).ppppresp = outcome(i).ppppstats.resp;
end
%
clear data data_left shaffle_p p

for i = 1:length(Sum)
    fprintf('Process Neurons %4.2f \n', i)
    unit.corr   = spike2eventRasteandPSTH_NP(outcome(i).pseucorLicks,outcome(i).dscorr,400,0,400);
    unit.error  = spike2eventRasteandPSTH_NP(outcome(i).pseuerrLicks,outcome(i).error,400,0,400);
    p(i) = Test_auROC_dISCRIMINATION(unit.corr.scmatrix,unit.error.scmatrix);
    % Using bootstrap to perform the statistic test
    com_data = [unit.corr .scmatrix; unit.error.scmatrix];
    for j = 1:1000 % 1000 iteration
        [data(:,j),idx_p] = datasample(com_data,length(unit.corr.scmatrix),'Replace',false);
        idx_pall          = 1:length(com_data);
        data_left(:,j)    = com_data(setdiff(idx_pall,idx_p));
        outcome(i).pppppLickstats.shaffle_p(j)      = Test_auROC_dISCRIMINATION(data(:,j),data_left(:,j));
    end
    if p(i)>=0 % one tail bootstrap
        outcome(i).pppppLickstats.p = length(find(outcome(i).pppppLickstats.shaffle_p>=p(i)))/1000;
    else
        outcome(i).pppppLickstats.p = length(find(outcome(i).pppppLickstats.shaffle_p<=p(i)))/1000;
    end
    
    if outcome(i).pppppLickstats.p>0.01
        outcome(i).pppppLickstats.resp = 0;
    else
        outcome(i).pppppLickstats.resp = 1;
    end
    clear data data_left shaffle_p 
end
% clear p
for i = 1:length(Sum)
    fprintf('Process Neurons %4.2f \n', i)
    unit.corr   = spike2eventRasteandPSTH_NP(Sum(i).timestampN,outcome(i).dscorr,400,0,400);
    unit.error  = spike2eventRasteandPSTH_NP(Sum(i).timestampN,outcome(i).error,400,0,400);
    p(i) = Test_auROC_dISCRIMINATION(unit.corr.scmatrix,unit.error.scmatrix);
    % Using bootstrap to perform the statistic test
    com_data = [unit.corr .scmatrix; unit.error.scmatrix];
    for j = 1:1000 % 1000 iteration
        [data(:,j),idx_p] = datasample(com_data,length(unit.corr.scmatrix),'Replace',false);
        idx_pall          = 1:length(com_data);
        data_left(:,j)    = com_data(setdiff(idx_pall,idx_p));
        outcome(i).pppppLickstats.shaffle_p(j)      = Test_auROC_dISCRIMINATION(data(:,j),data_left(:,j));
    end
    if p(i)>=0 % one tail bootstrap
        outcome(i).pppppstats.p = length(find(outcome(i).pppppLickstats.shaffle_p>=p(i)))/1000;
    else
        outcome(i).pppppstats.p = length(find(outcome(i).pppppLickstats.shaffle_p<=p(i)))/1000;
    end
    
    if outcome(i).pppppstats.p>0.01
        outcome(i).pppppstats.resp = 0;
    else
        outcome(i).pppppstats.resp = 1;
    end
    clear data data_left shaffle_p 
end
for i = 1:length(outcome)
    outcome(i).pppppresp = outcome(i).pppppstats.resp;
end
%%
save('outcome.mat','outcome')
%% Let's plot the curve of responsive outcome neurons
% p means 500 ms
% pp means 200 ms
% ppp means 100 ms
% pppp means 300 ms
% ppppp means 400 ms
field_name = {'ppp', 'pp', 'pppp', 'ppppp',''};
for i = 1 : 5
    quant_res(:,i) = [outcome.([field_name{i},'resp'])];
end
figure;
plot(sum(quant_res)./size(quant_res,1),'-o')
xlim([-1,6])
ylim([0,0.4])
% let's get the stats from error
for i = 1 : 5
    for j = 1: size(quant_res,1)
        if quant_res(j,i) ==0 
            quant_resLick(j,i) = 0;
        else
            quant_resLick(j,i) = outcome(j).([field_name{i},'Lickstats']).resp;
        end
    end
end
hold on
plot(sum(quant_resLick)./size(quant_resLick,1),'-o')
ylim([0,0.3])
%% let's get the duration of lateral lick for correct and error trials
mouse = unique({Sum.mouse});
for i = 1:length(mouse)
    temp = Sum(find(strcmp({Sum.mouse}, mouse{i})));
    date{i}  = unique({temp.date});
end
indx_forAna =[];
for i = 1:length(mouse)
    temp = date{i};
    for j = 1:length(temp)
        indx_mouse = find(strcmp({Sum.mouse}, mouse{i}));
        indx_date  = find(strcmp({Sum.date}, temp{j}));
        indx_common = intersect(indx_mouse,indx_date);
        indx_forAna = [indx_forAna, indx_common(1)];
    end
end

for i = 1:length(indx_forAna)
    licking = outcome(indx_forAna(i)).temp;
    indx_cor = find(cell2mat(licking(:,4)) ==1);
    indx_err = find(cell2mat(licking(:,4)) ==0);
    for j = 1:length(indx_cor)
        dura_cor(j) = licking{indx_cor(j),5}(end) -licking{indx_cor(j),5}(1);
    end
    for j = 1:length(indx_err)
        dura_err(j) = licking{indx_err(j),6}(end) -licking{indx_err(j),6}(1);
    end
    duration(i,1) = mean(dura_cor);
    duration(i,2) = mean(dura_err);
    clear dura_cor dura_err
end
%%
% plot the result and perform the statistic test
control = duration;
m1 = mean(control(:,1));
m2 = mean(control(:,2));
name = {[],'Correct','Error',[]}
sem1 = std(control(:,1))./sqrt(length(control(:,1)));
sem2 = std(control(:,2))./sqrt(length(control(:,2)));
figure;
bar(1,m1,'EdgeColor',[0 0 0],'FaceColor',[0.5,0.5,1],'LineWidth',1);hold on
bar(2,m2,'EdgeColor',[0 0 0],'FaceColor',[1,0.5,0.25],'LineWidth',1);hold on
set(gca,'XTick',0:3)
set(gca,'xticklabel',name)
ylabel('Duration (s)')
hold on;
e=errorbar([m1, m2], [sem1, sem2],'LineWidth',1.5);
e.Color ='k';
e.LineStyle = 'none';
ylim([0,3])
% perform the t-test
[h,p,~,stats] = ttest2(control(:,1),control(:,2));
%% let's identify which one is outcome
for i = 1:size(quant_res,1)
    for j = 1:size(quant_res,2)
        if quant_res(i,j) ==1 & quant_resLick(i,j) ==0
            outcome_resp(i,j) =1;
        else
             outcome_resp(i,j) =0;
            
        end
    end
end

length(find(sum(outcome_resp,2)>=1))
%% summarize
for i = 1:length(Sum)
    summary(i).planning = Sum(i).auROC.PlanRcorr_vsLcorr.stats.respFlag;
    summary(i).action   = Sum(i).auROC.Rcorr_vsLcorr.stats.respFlag;
    summary(i).outcome  = a(i);
end
%% save
save('summary.mat','summary')