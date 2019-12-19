load('Sum2.mat')
load('PlanningTaste.mat')
for i = 1:length(Sum)
    Idx(i) = Sum(i).auROC.PlanRcorr_vsLcorr.stats.respFlag;
    Idx_dir(i) = Sum(i).auROC.PlanRcorr_vsLcorr.p;
end

idd = find(Idx ==1); % find index of preparatory neurons
idd_t = idd(bothnon); % find the index of taste-selective preparatory neurons
idd_nont = setdiff(idd,idd_t); % find the index of exclusive preparatory neurons
Data_t_prep = Sum(idd_t);
Data_nont_prep = Sum(idd_nont);
save('data_prep.mat','Data_t_prep','Data_nont_prep')
%% reorganize data
clear
load('data_prep.mat')
Sum = Data_t_prep;
Sum = Data_nont_prep;
Sum = [Data_t_prep,Data_nont_prep];
%%
reorganize_data(Sum)
%
function reorganize_data(Sum)
for k = 1:length(Sum)
    spike = Sum(k).PSTH.RightRew.Corr.Spike;
    event.RCorr = Sum(k).event.tsRCorr.FLickRSpou;
%     event.LCorr = Sum(k).event.tsLCorr.FLickLSpou;
    event.RErr = Sum(k).event.tsRErr.FLickRSpou;
%     event.LErr  = Sum(k).event.tsLErr.FLickRSpou;


%     event.RCorr = Sum(k).event.tsRCorr.FLickTaste;
%     event.LCorr = Sum(k).event.tsLCorr.FLickTaste;
%     event.RErr = Sum(k).event.tsRErr.FLickTaste;
%     event.LErr  = Sum(k).event.tsLErr.FLickTaste;


%     event_sum = [event.RCorr; event.LCorr; event.RErr; event.LErr];
        event_sum = [event.RCorr;  event.RErr];
%             event_sum = [event.LCorr;  event.LErr];
    
        event_info= [ones(size(event.RCorr)); 2*ones(size(event.RErr))];
%         event_info= [ones(size(event.LCorr)); 2*ones(size(event.LErr))];
%     event_info= [ones(size(event.RCorr)); 2*ones(size(event.LCorr));3*ones(size(event.RErr)); 4*ones(size(event.LErr))];
    
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
%                 raster_labels.stimulus_ID{i} = 'LCorr';

%             case 3
%                 raster_labels.stimulus_ID{i} = 'RErr';
%             case 4
%                 raster_labels.stimulus_ID{i} = 'LErr';
        end
    end
    raster_site_info.alignment_event_time = 4501;
%     raster_site_info.alignment_event_time = 1501;

    clearvars -except raster_data raster_labels raster_site_info Sum k
    save(['Decision_neurodecoding_',num2str(k),'.mat'],'raster_data','raster_labels','raster_site_info')
    fprintf('Finish processing neuron # %0.1f\n', k)
    clearvars -except Sum k
end
end