function res=decoder_multi_bag_x_selectivity(spikes_train,spikes_test,options)
%
% Multi-class decoder
% Use bagging (bootstrap aggregating) to decrease variance of
% classification (Breitman, 1996).
% Use permutation test (trial-shuffled classification) to assess
% significance of performance
% Compute decoding accuracy for each stimulus, by summing euclidean
% distance over all time bins.
% 
% INPUT train_spikes,test_spikes =structure with fields=conditions (e.g. 'S','C','N','Q'). 
%                   Each field is a structure array of dim (ntrials,nunits) with field .spk
%                   containing array of spike times
%                   train and test spikes must contain the same number of
%                   conditions, train=training set,test=testing set
%       options=structure with fields:
%               .binsize=size of decoding bin (s) [for euclidean distance decoder only]
%               .window =size of decoding window
%               .xval   =k; k is the number of xval hold-outs (e.g.
%                        k=1 for 1 hold-res trial).
%               .mintrials (optional): specify smallest number of trials 
%                       (if field is absent, pick min across conditions.
%               .shuffled (optional): if 1, run also in shuffled mode; if 0
%                       or absent, no shuffled.
%
% OUPUT results=structure with fields 'data' and 'shuffled', 'pvalue' and 'pvalue_stimulus:
%               .data.accuracy=average accuracy in each bin
%               .data.confusion=confusion matrix, gives decoding accuracy with
%                           1st=training event; 2nd=testing event; 3rd=bin;
%               .pvalue=fraction of average shuffled accuracies>average empirical
%                               accuracy;
%               .pvalue_stimulus=fraction of shuffled accuracies>empirical
%                               accuracy for each stimulus (THIS IS THE IMPORTANT ONE);
%               .shuffled.ave_accuracy=shuffled accuracy averaged across stimuli
%               .shuffled.ave_accuracy_ci=95% CI of averaged shuffled distribution
%               .shuffled.accuracy_ci_perstimulus=95% CI of shuffled distribution for each taste
%               .shuffled.accuracy_boot_perstimulus=distribution of
%               shuffled accuracies by stimulus (used to compute
%               accuracy_ci_perstimulus and pvalue_stimulus) 
%               .shuffled.confusion=shuffled confusion matrix
%
% Luca Mazzucato August 2017

% DEBUG
% options=struct('binsize',0.2,'window',[-1 2],'xval',2);
% options.mintrials=9;
% spikes_train=spike;
% spikes_test=spike;
% options=opt;

% spikes_test=spikes_train;
trigger.train=fieldnames(spikes_train);
trigger.test=fieldnames(spikes_test);
numev=numel(trigger.train);
if numel(trigger.test)~=numev
    fprintf('error: # of conditions in training, testing sets are different!');
    return;
end

ntrials=struct('train',zeros(1,numev),'test',zeros(1,numev));
for E=1:numev
    ntrials.train(E)=size(spikes_train.(trigger.train{E}),1);
    ntrials.test(E)=size(spikes_test.(trigger.test{E}),1);
end
nunits=size(spikes_train.(trigger.train{E}),2);

% default options
type='pop';
binsize=0.1;
windows=[-1 1];
xval=2; % hold-res in cross-validation
nbag=1;
nboot=100;
% input
v2struct(options);
mintrials=min([ntrials.train ntrials.test]);
if any(strcmp(fieldnames(options),'mintrials'))
    mintrials=min(mintrials,options.mintrials);
end
nx=nchoosek(mintrials,xval); % number of xval runs
MODES={'data'};
% OUTPUT
res=struct('data',struct(),'nrun',nx,'ntrain',mintrials,'ntest',xval,'nbag',nbag,...
    'xbins',[]);    
%
if any(strcmp(fieldnames(options),'shuffled'))
    if shuffled==1
        MODES{2}='shuffled';
        res.shuffled=struct();
        res.pvalue=[];
        res.pvalue_stimulus=[];
    end
end



% BINS - trick to create sliding window
bin0=windows(1):binsize:windows(2);
% no moving window
slidewindow=bin0;
xbins=slidewindow(2:end);
[nsteps,nbins0]=size(xbins);
nbins=numel(xbins);
indbins=1:nbins;
res.xbins=xbins;

spkcnt=struct('train',zeros(mintrials,nbins,nunits,numev),'test',zeros(mintrials,nbins,nunits,numev));
for E=1:numev
    for st=1:nsteps
        % train
        temp=Spikes2Bins(spikes_train.(trigger.train{E}),slidewindow(st,:));
        a=randperm(ntrials.train(E));
%         a=1:ntrials.train(E);
        % only keep first mintrials trials, throw away the rest
        spkcnt.train(1:mintrials,(st-1)*nbins0+1:st*nbins0,1:nunits,E)=temp(a(1:mintrials),1:nbins0,1:nunits);
        % test
        temp=Spikes2Bins(spikes_test.(trigger.test{E}),slidewindow(st,:));
        if ~isequal(trigger.train,trigger.test)
            a=randperm(ntrials.test(E));
            spkcnt.test(1:mintrials,(st-1)*nbins0+1:st*nbins0,1:nunits,E)=temp(a(1:mintrials),1:nbins0,1:nunits);
        else % if train and test set are the same, do not reshuffle trials
            spkcnt.test=spkcnt.train;
        end
    end
end
spkcnt.train=spkcnt.train(:,indbins,:,:);
spkcnt.test=spkcnt.test(:,indbins,:,:);


%----------
% DECODING
%----------
% BAGGING
% bagging procedure for k-cross-validation
% 1) split the dataset into training set L and test set T (with k trials) for one step of
% the k-cross-validation procedure
% 2) for each training set, create B L^b bootstrapped sets b=1,...,B by
% sampling with replacement from L, and construct the classifier phi^b in
% each set (eg, the PSTH)
% 3) classify each x in the test set T with all classifiers phi^1,...phi^B,
% obtaining B different votes. The most frequent vote among the B votes is the bagged
% classification x.
% 4) repeat for each cross-validation step (L,T).
% 
% STATISTICAL SIGNIFICANCE
% Permutation test
% Repeat nboot times. For each run, shuffle all trials across classes and
% compute accuracy. The p-value of the data is the probability that the data
% exceeds the shuffled accuracy cdf.

psthfun=@(A)sum(A,1)/(mintrials-xval);
%------------
% XVAL TRIALS
%------------
% xval training and test sets, compare each pair
% trainlog and test set are structures whose fields trigger are nx *
% mintrials arrays: each row is a run, whose columns have 0s and 1s for
% training and testing trials in that run
indx=1:nx;
indxtrain=repmat(indx',1,mintrials-xval);
%
% pick mintrials from larger set for each bootstrap run
% % % % % databoot=repmat(struct('trainind',struct(),'testind',struct(),...
% % % % %     'spkcnt',struct('train',struct('shuffled',zeros(mintrials,nbins,nunits,numev),'data',[]),...
% % % % %                     'test',struct('shuffled',zeros(mintrials,nbins,nunits,numev),'data',[])),...
% % % % %     'template',zeros(nx,nbins,nunits,numev,nbag),...
% % % % %     'confusion',zeros(numev,numev,nbins)),1,nboot); % indices for training and testing for each bootstrap run
databoot=repmat(struct('trainind',struct(),'testind',struct(),...
    'spkcnt',struct('train',struct('shuffled',zeros(mintrials,nbins,nunits,numev),'data',[]),...
                    'test',struct('shuffled',zeros(mintrials,nbins,nunits,numev),'data',[])),...
    'template',zeros(nx,nbins,nunits,numev,nbag),...
    'confusion',zeros(numev,numev)),1,nboot); % indices for training and testing for each bootstrap run
% for shuffled, use boot=1:nboot, for data use boot=1; 
% for boot=1:nboot
parfor boot=1:nboot
    for E=1:numev
        % TRAINING AND TESTING INDICES for each bootstrap run
        % training
        indxval_all=nchoosek(1:mintrials,mintrials-xval); % all ways to pick mintrials-xval res of ntrials(E) trials
        a=randperm(nchoosek(mintrials,mintrials-xval)); % keep only nx random runs
        % train
        indxval_train=indxval_all(a(1:nx),1:mintrials-xval); % rows=runs; cols=indices of training trials in each run
%         indxval_train=indxval_all((1:nx),1:mintrials-xval); % rows=runs; cols=indices of training trials in each run
        databoot(boot).trainind.(trigger.train{E})=indxval_train;
        indonetrain=sub2ind([nx,mintrials],indxtrain,indxval_train); % linear index
%         trainlog.(trigger{E})=zeros(nx,ntrials(E));
%         trainlog.(trigger{E})(indonetrain)=1; % rows=runs; columns: 0s and 1s for training trials
        % test
        indxval_test=ones(nx,mintrials);
        indxval_test(indonetrain)=0; % all leftover trials from removing training set
        [I,J]=find(indxval_test);
        [~,indsort] = sort(I);
        JJ=J(indsort);
        indxval_test=reshape(JJ,[xval,nx])'; % all leftover trials (same format as indxval_train)
        indxval_test=indxval_test(:,1:xval); % keep only first xval leftover trials in each run
    %     indonetest=sub2ind([size(indxval_train,1),ntrials(E)],indxtest,indxval_test);
    %     testset.(trigger{E})=zeros(nx,ntrials(E));
    %     testset.(trigger{E})(indonetest)=1; % rows=run; columns: 0s and 1s for test trials
        databoot(boot).testind.(trigger.test{E})=indxval_test; % rows=run; columns: 0s and 1s for test trials
        % TRIALS FOR EACH BOOTSTRAP RUN sampled with replacement
%         indtrials=datasample(1:ntrials(E),mintrials);
        % nboot shuffled spike counts
        %
        % TRAIN
        spktemp=zeros(mintrials*numev,nbins,nunits);
        for E=1:numev
            spktemp((E-1)*mintrials+1:E*mintrials,1:nbins,1:nunits)=...
                spkcnt.train(1:mintrials,1:nbins,1:nunits,E);
        end
        a=randperm(numev*mintrials);
        spktemp=spktemp(a,:,:);
        for E=1:numev
            databoot(boot).spkcnt.train.shuffled(1:mintrials,1:nbins,1:nunits,E)=...
                spktemp((E-1)*mintrials+1:E*mintrials,1:nbins,1:nunits);
        end
        % TEST
        spktemp=zeros(mintrials*numev,nbins,nunits);
        for E=1:numev
            spktemp((E-1)*mintrials+1:E*mintrials,1:nbins,1:nunits)=...
                spkcnt.test(1:mintrials,1:nbins,1:nunits,E);
        end
        a=randperm(numev*mintrials);
        spktemp=spktemp(a,:,:);
        for E=1:numev
            databoot(boot).spkcnt.test.shuffled(1:mintrials,1:nbins,1:nunits,E)=...
                spktemp((E-1)*mintrials+1:E*mintrials,1:nbins,1:nunits);
        end
     end
end
% use databoot(1) indices for 'data', and all nboot indices the 'shuffled' 
for m=1:numel(MODES)
    % initialize variable for bootstrap
    switch MODES{m}
        case 'data'
            NBOOT=1;
            databoot(1).spkcnt.train.data=spkcnt.train;
            databoot(1).spkcnt.test.data=spkcnt.test;
        case 'shuffled'
            NBOOT=nboot;
    end
    confusion_boot=zeros(numev,numev,NBOOT);
    parfor boot=1:NBOOT
%     for boot=1:NBOOT;
        % initialize training and testing sets for XVAL
        % initialize each bootstrap run
        %--------------------
        % EUCLIDEAN DISTANCE
        %--------------------
        % bagging: for row in trainind, run 10 bootstrap templates'
        for E=1:numev
            for bi=1:nbins
                for u=1:nunits
                    temp=squeeze(databoot(boot).spkcnt.train.(MODES{m})(:,bi,u,E)); % spike count in all trials for current bin, unit, event
%                     BOOTSTAT=bootstrp(nbag,@(x)psthfun(x),temp(databoot(boot).trainind.(trigger.train{E}))'); % dim=(nbag,nx); psth in each bag and xval run
                    BOOTSTAT=psthfun(temp(databoot(boot).trainind.(trigger.train{E}))'); % dim=(nbag,nx); psth in each bag and xval run
                    databoot(boot).template(1:nx,bi,u,E,1:nbag)=BOOTSTAT';
                end
            end
        end
        % 
        % test trials
        for Etest=1:numev
            % for each test trial (1st) get the eucl. dist. from the numev templates (4th) for
            % each bin (2nd) and unit (3rd)
            dist=zeros(xval*nx,numev,nbag);
            for x=1:xval
                % test dim (nx,bins,unit);
                indtest=databoot(boot).testind.(trigger.test{Etest})(:,x);
                test=squeeze(databoot(boot).spkcnt.test.(MODES{m})(indtest,:,:,Etest));
                stest=size(test);
                temp_test=repmat(test,[ones(1,numel(stest)),numev, nbag]);
                if nunits==1
                    % eucl dist without sum over neurons: (nx,numev,nbag)
                    temp_template=squeeze(databoot(boot).template); %(nx,nbins,numev,nbag)
                    dist((x-1)*nx+1:x*nx,1:numev,1:nbag)=...
                        sum((temp_test-temp_template).^2,2); % (nx,numev,nbag) sum over bins!
% % % % %                 else
% % % % %                     % eucl dist summed over all neurons (nx,nbins,numev,nbag)
% % % % %                     dist((x-1)*nx+1:x*nx,1:nbins,1:numev,1:nbag)=...
% % % % %                         squeeze(sum((temp_test-databoot(boot).template).^2,3));
                end
            end
            % bootstrap votes
            % dist: [xval runs, numev, nbag]
            [Y,indbags]=min(dist,[],2); % votes for each bag
            for k=2:numev
                for i_bag=1:nbag
                    tie=squeeze(sum(dist(:,:,i_bag)==repmat(Y(:,:,i_bag),1,numev,1),2)==k); % (xval,1)
                    I=find(tie); % how many ties
                    if ~isempty(I)
                        for i_tie=1:numel(I)
                            indtie=squeeze(dist(I(i_tie),:,i_bag)==Y(I(i_tie),:,i_bag));
                            itemp=indtie.*rand(size(indtie));
                            [~, bind]=max(itemp);
                            indbags(I(i_tie),1,i_bag)=bind;
                        end
                    end
                end
            end
            % for each test trial, get most frequent vote across
            % bags; when there are ties in the mode, adjudicate to first label
            % (Breitman, 1996)
            ind = mode(indbags,3); % (nx*xval,nbins) array of votes
            for Etrain=1:numev
                databoot(boot).confusion(Etest,Etrain,1:nbins)=sum(ind==Etrain*ones(size(ind)),1)/(nx*xval);
            end
        end    
    end
    % reformat confusion matrix 
    for boot=1:NBOOT
        confusion_boot(1:numev,1:numev,boot)=...
            databoot(boot).confusion(1:numev,1:numev);
    end
    % OUTPUT
    % pick diag of confusion matrix
    tempacc_stimulus=NaN(numev,NBOOT);
    for E=1:numev
        tempacc_stimulus(E,1:NBOOT)=squeeze(confusion_boot(E,E,1:NBOOT))';
    end
    tempacc=squeeze(mean(tempacc_stimulus,1));
    acc.(MODES{m})=tempacc; % NBOOT x nbins
    acc_stimulus.(MODES{m})=tempacc_stimulus; % NBOOT 
    confusion=mean(confusion_boot,3); % numev x numev 
    % output variables
    if strcmp(MODES{m},'shuffled')
        accuracy_ci_stimulus=NaN(numev,2);
        accuracy=squeeze(mean(acc.(MODES{m}),1));
        accuracy_ci=prctile(acc.(MODES{m}),[2.5,97.5],1);
        for E=1:numev
            temp_prc=prctile(tempacc_stimulus(E,1:NBOOT),[2.5,97.5],2);
            accuracy_ci_stimulus(E,1:2)=temp_prc(1,1:2); %nd changed from "temp_prc(1:2) to 3d matrix
        end
        res.(MODES{m})=struct('ave_accuracy',accuracy,'ave_accuracy_ci',accuracy_ci,'accuracy_ci_perstimulus',accuracy_ci_stimulus,...
            'confusion',confusion,'nboot',NBOOT,'ave_accuracy_boot',tempacc,'accuracy_boot_perstimulus',tempacc_stimulus);    
    else
        accuracy=acc.(MODES{m})';
        res.(MODES{m})=struct('accuracy',accuracy,...
            'confusion',confusion);    
    end
end

% IF SHUFFLE MODE TEST SIGNIFICANCE
if shuffled
    pvalue=NaN(1);
    pvalue_stimulus=NaN(numev,1);
    % pick diag of confusion matrix
    empirical_cdf=acc.shuffled(:);
    pvalue=sum(empirical_cdf>acc.data)/numel(empirical_cdf);
    for E=1:numev
        empirical_cdf=squeeze(acc_stimulus.shuffled(E,:));
        pvalue_stimulus(E)=sum(empirical_cdf>res.data.confusion(E,E))/numel(empirical_cdf);
    end
    res.pvalue=pvalue;
    res.pvalue_stimulus=pvalue_stimulus;
end

