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