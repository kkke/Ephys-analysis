%% load data
load('Sum2.mat')
for i = 1:length(Sum) % reorganize the data
    tastepsth(i).M = Sum(i).PSTH.RightRew.Maltose;
    tastepsth(i).O = Sum(i).PSTH.RightRew.Octacetate;
    tastepsth(i).S = Sum(i).PSTH.LeftRew.Sucrose;
    tastepsth(i).Q = Sum(i).PSTH.LeftRew.Quinine;
end
% reorganize the data
taste = {'S', 'M', 'Q','O'};
for i = 1:length(Sum) % recalculate the psth with 50 ms as bin size, and only take from 0 to 500 ms
    for j = 1:length(taste)
        tastepsth_bin25(i).(taste{j}) = spike2eventRasteandPSTH_NP(tastepsth(i).(taste{j}).Spike, tastepsth(i).(taste{j}).Event, 50,0,2500);
    end
    fprintf('Finish processing neuron # %0.f\n',i)
end
%% start single neuron decoding analysis
% Use server to perform it
% for i = 1:length(Sum)
%     spikes.S = tastepsth_bin25(i).S.spikeraster';
%     spikes.M = tastepsth_bin25(i).M.spikeraster';
%     spikes.Q = tastepsth_bin25(i).Q.spikeraster';
%     spikes.O = tastepsth_bin25(i).O.spikeraster';
%     %%
%     trig=fieldnames(spikes);
%     nev=numel(trig);
%     %%
%     % decoding parameters
%     kind='units';
%     xval=2; % number of leave-out trials for each x-validation run
%     binsize=0.05; % size of decoding window
%     % windows=[0 2.5]; % total trial interval to be decoded aligned at t=0.
%     windows=[0 2.5]; % total trial interval to be decoded aligned at t=0.
%     nboot=1000; % number of bootstrap runs for each shuffled stimulus
% %     filesave_dec='decoding'; % prefix of files to be saved
%     %%
%     [results(i), options(i)] = fun_decode_units_selectivity(spikes,binsize,windows,xval,nboot,trig);
%     fprintf('Finish decoding neuron # %0.1f \n', i)
% %     fprintf
% end
%% check for each individual neurons

i = 9
trig = {'S','M','Q','O'};
color = {'r','b','g','k'};
fun_decode_units_selectivity_plot(results(i),trig')
figure
ok =1;
for j = 1:length(trig)
    spiketime = tastepsth_bin25(i).(trig{j}).spikeraster;
    subplot(2,1,1)
    for k = 1:length(spiketime)
        if isempty(spiketime(k).times)
        else
            scatter(spiketime(k).times, ok*ones(size(spiketime(k).times)),[color{j},'.'])
            hold on
            ok = ok+1;
        end
    end
    subplot(2,1,2)
    plot(tastepsth_bin25(i).(trig{j}).timepoint, tastepsth_bin25(i).(trig{j}).FR_avg,[color{j},'-'])
    hold on
end

%% extract neurons with significant taste selectivity
for i = 1:length(results)
    results(i).selectivity = find(results(i).all.pvalue_stimulus< 0.05/4);
    if isempty(results(i).selectivity)
        results(i).sele_flag = 0;
    elseif results(i).all.pvalue>=0.05
        results(i).sele_flag = 0;
    else
        results(i).sele_flag = 1;
    end
end





