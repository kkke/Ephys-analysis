% OTHER PARAMETERS FOR DECODING
nbag=1; % for bagging, recommended if the number of trials is < 10
shuffled=1; % 1 to include shuffle analysis
fieldNames={'binsize','window','xval','nboot','nbag','shuffled','fieldNames'};
options=v2struct(fieldNames);
%

% initialize empty variables
results.all=struct('data',[],...
    'nrun',[],'ntrain',[],'ntest',xval,'nbag',nbag,'xbins',[],'shuffled',[],'pvalue',[],'pvalue_stimulus',[]);

% find smallest number of trials across events
ntrials=zeros(1,numel(trig));
for E=1:numel(trig)
    [ntrials(E), nunits]=size(spikes.(trig{E}));
end
options.mintrials=min(ntrials); % pick smallest number of trials across conditions for decoding\

tic
spikes_train=spikes; spikes_test=spikes;
results.all=decoder_multi_bag_x_selectivity(spikes_train,spikes_test,options);
toc
filesave=[filesave_dec '_selectivity.mat'];
save(filesave,'results','options');
fprintf('Data saved to %s\n',filesave);
