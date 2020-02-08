%% reorganize data into raster format for both correc and error trials
% Sum = data_taste;
for k = 1:length(Sum)
    fprintf('Processing # %0.01f \n',k)
    spike = Sum(k).PSTH.RightRew.Corr.Spike;
    event.S = Sum(k).event.tsLCorr.TasteID.S;
    event.Q = Sum(k).event.tsLCorr.TasteID.Q;
    first_SQ = Sum(k).event.tsLCorr.TasteDel(1,1);
    if abs(first_SQ-event.S(1))<abs(first_SQ-event.Q(1))
        id_s = Sum(k).event.tsLCorr.TasteDel(1,2);
        id_q = setdiff(unique(Sum(k).event.tsLCorr.TasteDel(:,2)),id_s);
    else
        id_q = Sum(k).event.tsLCorr.TasteDel(1,2);
        id_s = setdiff(unique(Sum(k).event.tsLCorr.TasteDel(:,2)),id_q);
    end
    event.S_r = Sum(k).event.tsLErr.FLickTaste(find(Sum(k).event.tsLErr.TasteDel(:,2) == id_s));
    event.Q_r = Sum(k).event.tsLErr.FLickTaste(find(Sum(k).event.tsLErr.TasteDel(:,2) == id_q));
    
    event.M = Sum(k).event.tsRCorr.TasteID.M;
    event.SO =Sum(k).event.tsRCorr.TasteID.O;
    first_MO = Sum(k).event.tsRCorr.TasteDel(1,1);
    if abs(first_MO-event.M(1))<abs(first_MO-event.SO(1))
        id_m = Sum(k).event.tsRCorr.TasteDel(1,2);
        id_so = setdiff(unique(Sum(k).event.tsRCorr.TasteDel(:,2)),id_m);
    else
        id_so = Sum(k).event.tsRCorr.TasteDel(1,2);
        id_m = setdiff(unique(Sum(k).event.tsRCorr.TasteDel(:,2)),id_so);
    end
    event.M_r = Sum(k).event.tsRErr.FLickTaste(find(Sum(k).event.tsRErr.TasteDel(:,2) == id_m));
    event.SO_r = Sum(k).event.tsRErr.FLickTaste(find(Sum(k).event.tsRErr.TasteDel(:,2) == id_so)); 
    
    
    event_sum = [event.S; event.Q; event.M; event.SO; event.S_r; event.Q_r; event.M_r; event.SO_r];
    event_info= [ones(size(event.S)); 2*ones(size(event.Q)); 3*ones(size(event.M)); 4*ones(size(event.SO));...
        5*ones(size(event.S_r)); 6*ones(size(event.Q_r)); 7*ones(size(event.M_r)); 8*ones(size(event.SO_r))];
    event_sum(:,2) = event_info;
    [~, I] = sort(event_sum(:,1));
    for i = 1:length(I)
        event_order(i,:) = event_sum(I(i),:);
    end
    unit = spike2eventRasteandPSTH_NP(spike,event_order,1,-1500,2500); % use 1 ms to construct the raster format
    raster_data = unit.scmatrix;
    for i = 1:length(event_order)
        switch event_order(i,2)
            case 1
                raster_labels.stimulus_ID{i} = 'S';
            case 2
                raster_labels.stimulus_ID{i} = 'Q';
            case 3
                raster_labels.stimulus_ID{i} = 'M';
            case 4
                raster_labels.stimulus_ID{i} = 'SO';
            case 5
                raster_labels.stimulus_ID{i} = 'S_r';
            case 6
                raster_labels.stimulus_ID{i} = 'Q_r';
            case 7
                raster_labels.stimulus_ID{i} = 'M_r';
            case 8
                raster_labels.stimulus_ID{i} = 'SO_r';
        end
    end
    raster_site_info.alignment_event_time = 1501;
    clearvars -except raster_data raster_labels raster_site_info Sum k
    save(['Decision_neurodecoding_',num2str(k),'.mat'],'raster_data','raster_labels','raster_site_info')
    clearvars -except Sum k
end
%%

%% reorganize data into raster format for both correc and error trials
% Sum = data_taste;
for k = 1:length(Sum)
    fprintf('Processing # %0.01f \n',k)
    spike = Sum(k).PSTH.RightRew.Corr.Spike;
    event.S = Sum(k).event.tsLCorr.TasteID.S;
    event.Q = Sum(k).event.tsLCorr.TasteID.Q;
    first_SQ = Sum(k).event.tsLCorr.TasteDel(1,1);
    if abs(first_SQ-event.S(1))<abs(first_SQ-event.Q(1))
        id_s = Sum(k).event.tsLCorr.TasteDel(1,2);
        id_q = setdiff(unique(Sum(k).event.tsLCorr.TasteDel(:,2)),id_s);
    else
        id_q = Sum(k).event.tsLCorr.TasteDel(1,2);
        id_s = setdiff(unique(Sum(k).event.tsLCorr.TasteDel(:,2)),id_q);
    end
    event.S_r = Sum(k).event.tsLErr.FLickTaste(find(Sum(k).event.tsLErr.TasteDel(:,2) == id_s));
    event.Q_r = Sum(k).event.tsLErr.FLickTaste(find(Sum(k).event.tsLErr.TasteDel(:,2) == id_q));
    
    event.M = Sum(k).event.tsRCorr.TasteID.M;
    event.SO =Sum(k).event.tsRCorr.TasteID.O;
    first_MO = Sum(k).event.tsRCorr.TasteDel(1,1);
    if abs(first_MO-event.M(1))<abs(first_MO-event.SO(1))
        id_m = Sum(k).event.tsRCorr.TasteDel(1,2);
        id_so = setdiff(unique(Sum(k).event.tsRCorr.TasteDel(:,2)),id_m);
    else
        id_so = Sum(k).event.tsRCorr.TasteDel(1,2);
        id_m = setdiff(unique(Sum(k).event.tsRCorr.TasteDel(:,2)),id_so);
    end
    event.M_r = Sum(k).event.tsRErr.FLickTaste(find(Sum(k).event.tsRErr.TasteDel(:,2) == id_m));
    event.SO_r = Sum(k).event.tsRErr.FLickTaste(find(Sum(k).event.tsRErr.TasteDel(:,2) == id_so)); 
    
    event.S_d = randsample(event.S,length(event.S_r));
    if length(event.Q)<=length(event.Q_r)
        event.Q_d = event.Q
    else
        event.Q_d = randsample(event.Q,length(event.Q_r));
    end
    event.M_d = randsample(event.M,length(event.M_r));
    event.SO_d = randsample(event.S,length(event.SO_r));
    
    event_sum = [event.S_d; event.Q_d; event.M_d; event.SO_d];

 
    
    event_info= [ones(size(event.S_d)); 2*ones(size(event.Q_d)); 3*ones(size(event.M_d)); 4*ones(size(event.SO_d))];
    event_sum(:,2) = event_info;
    [~, I] = sort(event_sum(:,1));
    for i = 1:length(I)
        event_order(i,:) = event_sum(I(i),:);
    end
    unit = spike2eventRasteandPSTH_NP(spike,event_order,1,-1500,2500); % use 1 ms to construct the raster format
    raster_data = unit.scmatrix;
    for i = 1:length(event_order)
        switch event_order(i,2)
            case 1
                raster_labels.stimulus_ID{i} = 'S_d';
            case 2
                raster_labels.stimulus_ID{i} = 'Q_d';
            case 3
                raster_labels.stimulus_ID{i} = 'M_d';
            case 4
                raster_labels.stimulus_ID{i} = 'SO_d';
        end
    end
    raster_site_info.alignment_event_time = 1501;
    clearvars -except raster_data raster_labels raster_site_info Sum k
    save(['Decision_neurodecoding_',num2str(k),'.mat'],'raster_data','raster_labels','raster_site_info')
    clearvars -except Sum k
end
%%
%% reorganize data into raster format for both correc and error trials
% Sum = data_taste;
for k = 1:length(Sum)
    fprintf('Processing # %0.01f \n',k)
    spike = Sum(k).PSTH.RightRew.Corr.Spike;
    event.S = Sum(k).event.tsLCorr.TasteID.S;
    event.Q = Sum(k).event.tsLCorr.TasteID.Q;
    first_SQ = Sum(k).event.tsLCorr.TasteDel(1,1);
    if abs(first_SQ-event.S(1))<abs(first_SQ-event.Q(1))
        id_s = Sum(k).event.tsLCorr.TasteDel(1,2);
        id_q = setdiff(unique(Sum(k).event.tsLCorr.TasteDel(:,2)),id_s);
    else
        id_q = Sum(k).event.tsLCorr.TasteDel(1,2);
        id_s = setdiff(unique(Sum(k).event.tsLCorr.TasteDel(:,2)),id_q);
    end
    event.S_r = Sum(k).event.tsLErr.FLickTaste(find(Sum(k).event.tsLErr.TasteDel(:,2) == id_s));
    event.Q_r = Sum(k).event.tsLErr.FLickTaste(find(Sum(k).event.tsLErr.TasteDel(:,2) == id_q));
    
    event.M = Sum(k).event.tsRCorr.TasteID.M;
    event.SO =Sum(k).event.tsRCorr.TasteID.O;
    first_MO = Sum(k).event.tsRCorr.TasteDel(1,1);
    if abs(first_MO-event.M(1))<abs(first_MO-event.SO(1))
        id_m = Sum(k).event.tsRCorr.TasteDel(1,2);
        id_so = setdiff(unique(Sum(k).event.tsRCorr.TasteDel(:,2)),id_m);
    else
        id_so = Sum(k).event.tsRCorr.TasteDel(1,2);
        id_m = setdiff(unique(Sum(k).event.tsRCorr.TasteDel(:,2)),id_so);
    end
    event.M_r = Sum(k).event.tsRErr.FLickTaste(find(Sum(k).event.tsRErr.TasteDel(:,2) == id_m));
    event.SO_r = Sum(k).event.tsRErr.FLickTaste(find(Sum(k).event.tsRErr.TasteDel(:,2) == id_so)); 
    
    
    event_sum = [event.S; event.Q; event.M; event.SO; event.S_r; event.Q_r; event.M_r; event.SO_r];
    event_info= [ones(size(event.S)); 2*ones(size(event.Q)); 3*ones(size(event.M)); 4*ones(size(event.SO));...
        5*ones(size(event.S_r)); 6*ones(size(event.Q_r)); 7*ones(size(event.M_r)); 8*ones(size(event.SO_r))];
    event_sum(:,2) = event_info;
    [~, I] = sort(event_sum(:,1));
    for i = 1:length(I)
        event_order(i,:) = event_sum(I(i),:);
    end
    unit = spike2eventRasteandPSTH_NP(spike,event_order,1,-1500,2500); % use 1 ms to construct the raster format
    raster_data = unit.scmatrix;
    for i = 1:length(event_order)
        switch event_order(i,2)
            case 1
                raster_labels.stimulus_ID{i} = 'S';
            case 2
                raster_labels.stimulus_ID{i} = 'Q';
            case 3
                raster_labels.stimulus_ID{i} = 'M';
            case 4
                raster_labels.stimulus_ID{i} = 'SO';
            case 5
                raster_labels.stimulus_ID{i} = 'S_r';
            case 6
                raster_labels.stimulus_ID{i} = 'Q_r';
            case 7
                raster_labels.stimulus_ID{i} = 'M_r';
            case 8
                raster_labels.stimulus_ID{i} = 'SO_r';
        end
    end
    raster_site_info.alignment_event_time = 1501;
    clearvars -except raster_data raster_labels raster_site_info Sum k
    save(['Decision_neurodecoding_',num2str(k),'.mat'],'raster_data','raster_labels','raster_site_info')
    clearvars -except Sum k
end
%%

%% reorganize data into raster format for both correc and error trials? but just for Left-cued trials
% Sum = data_taste;
for k = 1:length(Sum)
    fprintf('Processing # %0.01f \n',k)
    spike = Sum(k).PSTH.RightRew.Corr.Spike;
    event.S = Sum(k).event.tsLCorr.TasteID.S;
    event.Q = Sum(k).event.tsLCorr.TasteID.Q;
    first_SQ = Sum(k).event.tsLCorr.TasteDel(1,1);
    if abs(first_SQ-event.S(1))<abs(first_SQ-event.Q(1))
        id_s = Sum(k).event.tsLCorr.TasteDel(1,2);
        id_q = setdiff(unique(Sum(k).event.tsLCorr.TasteDel(:,2)),id_s);
    else
        id_q = Sum(k).event.tsLCorr.TasteDel(1,2);
        id_s = setdiff(unique(Sum(k).event.tsLCorr.TasteDel(:,2)),id_q);
    end
    event.S_r = Sum(k).event.tsLErr.FLickTaste(find(Sum(k).event.tsLErr.TasteDel(:,2) == id_s));
    event.Q_r = Sum(k).event.tsLErr.FLickTaste(find(Sum(k).event.tsLErr.TasteDel(:,2) == id_q));
    
    event.M = Sum(k).event.tsRCorr.TasteID.M;
    event.SO =Sum(k).event.tsRCorr.TasteID.O;
    first_MO = Sum(k).event.tsRCorr.TasteDel(1,1);
    if abs(first_MO-event.M(1))<abs(first_MO-event.SO(1))
        id_m = Sum(k).event.tsRCorr.TasteDel(1,2);
        id_so = setdiff(unique(Sum(k).event.tsRCorr.TasteDel(:,2)),id_m);
    else
        id_so = Sum(k).event.tsRCorr.TasteDel(1,2);
        id_m = setdiff(unique(Sum(k).event.tsRCorr.TasteDel(:,2)),id_so);
    end
    event.M_r = Sum(k).event.tsRErr.FLickTaste(find(Sum(k).event.tsRErr.TasteDel(:,2) == id_m));
    event.SO_r = Sum(k).event.tsRErr.FLickTaste(find(Sum(k).event.tsRErr.TasteDel(:,2) == id_so)); 
    
    event.S_d = randsample(event.S,length(event.S_r));
    if length(event.Q)<=length(event.Q_r)
        event.Q_d = event.Q
    else
        event.Q_d = randsample(event.Q,length(event.Q_r));
    end
    event.M_d = randsample(event.M,length(event.M_r));
    event.SO_d = randsample(event.S,length(event.SO_r));
    
    %
    event_sum = [event.M_d; event.SO_d];
    event_info= [ones(size(event.M_d)); 2*ones(size(event.SO_d))];
    if length(event_sum)<2
        continue
    else
    end
    event_sum(:,2) = event_info;
    [~, I] = sort(event_sum(:,1));
    for i = 1:length(I)
        event_order(i,:) = event_sum(I(i),:);
    end
    unit = spike2eventRasteandPSTH_NP(spike,event_order,1,-1500,2500); % use 1 ms to construct the raster format
    raster_data = unit.scmatrix;
    for i = 1:length(event_order)
        switch event_order(i,2)
            case 1
                raster_labels.stimulus_ID{i} = 'M_d';
            case 2
                raster_labels.stimulus_ID{i} = 'SO_d';
            case 3
                raster_labels.stimulus_ID{i} = 'M_d';
            case 4
                raster_labels.stimulus_ID{i} = 'SO_d';
        end
    end
    raster_site_info.alignment_event_time = 1501;
    clearvars -except raster_data raster_labels raster_site_info Sum k
    save(['Decision_neurodecoding_',num2str(k),'.mat'],'raster_data','raster_labels','raster_site_info')
    clearvars -except Sum k
end
%% for delay align to the first lateral lick
for k = 1:length(Sum)
    fprintf('Processing # %0.01f \n',k)
    spike = Sum(k).PSTH.RightRew.Corr.Spike;
    event.S = Sum(k).event.tsLCorr.TasteID.S;
    event.Q = Sum(k).event.tsLCorr.TasteID.Q;
    first_SQ = Sum(k).event.tsLCorr.TasteDel(1,1);
    if abs(first_SQ-event.S(1))<abs(first_SQ-event.Q(1))
        id_s = Sum(k).event.tsLCorr.TasteDel(1,2);
        id_q = setdiff(unique(Sum(k).event.tsLCorr.TasteDel(:,2)),id_s);
    else
        id_q = Sum(k).event.tsLCorr.TasteDel(1,2);
        id_s = setdiff(unique(Sum(k).event.tsLCorr.TasteDel(:,2)),id_q);
    end
    event.S_r = Sum(k).event.tsLErr.FLickRSpou(find(Sum(k).event.tsLErr.TasteDel(:,2) == id_s));
    event.Q_r = Sum(k).event.tsLErr.FLickRSpou(find(Sum(k).event.tsLErr.TasteDel(:,2) == id_q));
    
    event.M = Sum(k).event.tsRCorr.TasteID.M;
    event.SO =Sum(k).event.tsRCorr.TasteID.O;
    first_MO = Sum(k).event.tsRCorr.TasteDel(1,1);
    if abs(first_MO-event.M(1))<abs(first_MO-event.SO(1))
        id_m = Sum(k).event.tsRCorr.TasteDel(1,2);
        id_so = setdiff(unique(Sum(k).event.tsRCorr.TasteDel(:,2)),id_m);
    else
        id_so = Sum(k).event.tsRCorr.TasteDel(1,2);
        id_m = setdiff(unique(Sum(k).event.tsRCorr.TasteDel(:,2)),id_so);
    end
    event.M_r = Sum(k).event.tsRErr.FLickRSpou(find(Sum(k).event.tsRErr.TasteDel(:,2) == id_m));
    event.SO_r = Sum(k).event.tsRErr.FLickRSpou(find(Sum(k).event.tsRErr.TasteDel(:,2) == id_so)); 
    
    event.S_d = randsample(event.S,length(event.S_r));
    if length(event.Q)<=length(event.Q_r)
        event.Q_d = event.Q
    else
        event.Q_d = randsample(event.Q,length(event.Q_r));
    end
    event.M_d = randsample(event.M,length(event.M_r));
    event.SO_d = randsample(event.S,length(event.SO_r));
    
    %
    event_sum = [event.M_r; event.SO_r];
    event_info= [ones(size(event.M_r)); 2*ones(size(event.SO_r))];
    if length(event_sum)<2
        continue
    else
    end
    event_sum(:,2) = event_info;
    [~, I] = sort(event_sum(:,1));
    for i = 1:length(I)
        event_order(i,:) = event_sum(I(i),:);
    end
    unit = spike2eventRasteandPSTH_NP(spike,event_order,1,-2000,1000); % use 1 ms to construct the raster format
    raster_data = unit.scmatrix;
    for i = 1:length(event_order)
        switch event_order(i,2)
            case 1
                raster_labels.stimulus_ID{i} = 'M_r';
            case 2
                raster_labels.stimulus_ID{i} = 'SO_r';
            case 3
                raster_labels.stimulus_ID{i} = 'M_d';
            case 4
                raster_labels.stimulus_ID{i} = 'SO_d';
        end
    end
    raster_site_info.alignment_event_time = 2001;
    clearvars -except raster_data raster_labels raster_site_info Sum k
    save(['Decision_neurodecoding_',num2str(k),'.mat'],'raster_data','raster_labels','raster_site_info')
    clearvars -except Sum k
end