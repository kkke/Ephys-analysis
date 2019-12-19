%% reorganize data into raster format
% Sum = data_taste;
for k = 1:length(Sum)
    spike = Sum(k).PSTH.RightRew.Corr.Spike;
    event.RCorr = Sum(k).event.tsRCorr.FLickRSpou;
%     event.LCorr = Sum(k).event.tsLCorr.FLickLSpou;
    event.RErr = Sum(k).event.tsRErr.FLickRSpou;
%     event.LErr  = Sum(k).event.tsLErr.FLickRSpou;
%     event_sum = [event.RCorr; event.LCorr; event.RErr; event.LErr];
    event_sum = [event.RCorr;  event.RErr];
%     event_sum = [event.LCorr;  event.LErr];
    
    event_info= [ones(size(event.RCorr)); 2*ones(size(event.RErr))];
    
%     event_info= [ones(size(event.LCorr)); 2*ones(size(event.LErr))];

    event_sum(:,2) = event_info;
    [~, I] = sort(event_sum(:,1));
    for i = 1:length(I)
        event_order(i,:) = event_sum(I(i),:);
    end
    unit = spike2eventRasteandPSTH_NP(spike,event_order,1,-4500,1500); % use 1 ms to construct the raster format
    raster_data = unit.scmatrix;
    for i = 1:length(event_order)
        switch event_order(i,2)
            case 1
                raster_labels.stimulus_ID{i} = 'RCorr';
%             case 1
%                 raster_labels.stimulus_ID{i} = 'LCorr';
            case 2
                raster_labels.stimulus_ID{i} = 'RErr';
%             case 2
%                 raster_labels.stimulus_ID{i} = 'LErr';
        end
    end
    raster_site_info.alignment_event_time = 4501;
    clearvars -except raster_data raster_labels raster_site_info Sum k
    save(['Decision_neurodecoding_',num2str(k),'.mat'],'raster_data','raster_labels','raster_site_info')
    fprintf('Finish processing neuron # %0.1f\n', k)
    clearvars -except Sum k
end
%%