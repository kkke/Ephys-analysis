%Use server to perform it
for i = 1:length(Sum)
    spikes.S = tastepsth_bin25(i).S.spikeraster';
    spikes.M = tastepsth_bin25(i).M.spikeraster';
    spikes.Q = tastepsth_bin25(i).Q.spikeraster';
    spikes.O = tastepsth_bin25(i).O.spikeraster';
    %%
    trig=fieldnames(spikes);
    nev=numel(trig);
    %%
    % decoding parameters
    kind='units';
    xval=2; % number of leave-out trials for each x-validation run
    binsize=0.05; % size of decoding window
    % windows=[0 2.5]; % total trial interval to be decoded aligned at t=0.
    f = {'bin0_05','bin05_10','bin10_15','bin15_20','bin20_25'};
    for j = 1:5
        windows=[j*0.5-0.5, j*0.5]; % total trial interval to be decoded aligned at t=0.
        nboot=1000; % number of bootstrap runs for each shuffled stimulus
        %     filesave_dec='decoding'; % prefix of files to be saved
        %%
        [results(i).(f{j}), ~] = fun_decode_units_selectivity(spikes,binsize,windows,xval,nboot,trig);
        fprintf('Finish decoding neuron # %0.1f \n', i)
    end
%     fprintf
end