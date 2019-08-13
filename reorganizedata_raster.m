%% reorganize data into raster format
Sum = data_taste;
for k = 1:length(Sum)
    spike = Sum(k).PSTH.RightRew.Corr.Spike;
    event.S = Sum(k).PSTH.LeftRew.Sucrose.Event;
    event.Q = Sum(k).PSTH.LeftRew.Quinine.Event;
    event.M = Sum(k).PSTH.RightRew.Maltose.Event;
    event.SO =Sum(k).PSTH.RightRew.Octacetate.Event;
    event_sum = [event.S; event.Q; event.M; event.SO];
    event_info= [ones(size(event.S)); 2*ones(size(event.Q)); 3*ones(size(event.M)); 4*ones(size(event.SO))];
    event_sum(:,2) = event_info;
    [~, I] = sort(event_sum(:,1));
    for i = 1:length(I)
        event_order(i,:) = event_sum(I(i),:);
    end
    unit = spike2eventRasteandPSTH_NP(spike,event_order,100,-1500,2500);
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
        end
    end
    raster_site_info.alignment_event_time = 1501;
    clearvars -except raster_data raster_labels raster_site_info Sum k
    save(['Decision_neurodecoding_',num2str(k),'.mat'],'raster_data','raster_labels','raster_site_info')
    clearvars -except Sum k
end
%%