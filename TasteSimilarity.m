load('Sum2.mat')
%%
% Taste Responsivness calculated with Ranksum (baseline -500ms Vs. Evoked
% 500ms); firing rate value from all trials and averaged across bins
for k = 1:size(Sum,2)
    tasteResp2(k,1) = ranksum(mean(Sum(k).PSTH.RightRew.Maltose.FRmatrix(:,5:10),2),...
        mean(Sum(k).PSTH.RightRew.Maltose.FRmatrix(:,11:16),2)); % ranksum test for each taste
    tasteResp2(k,2) = ranksum(mean(Sum(k).PSTH.RightRew.Octacetate.FRmatrix(:,5:10),2),...
        mean(Sum(k).PSTH.RightRew.Octacetate.FRmatrix(:,11:16),2));
    tasteResp2(k,3) = ranksum(mean(Sum(k).PSTH.LeftRew.Sucrose.FRmatrix(:,5:10),2),...
        mean(Sum(k).PSTH.LeftRew.Sucrose.FRmatrix(:,11:16),2));
    tasteResp2(k,4) = ranksum(mean(Sum(k).PSTH.LeftRew.Quinine.FRmatrix(:,5:10),2),...
        mean(Sum(k).PSTH.LeftRew.Quinine.FRmatrix(:,11:16),2));
    for t=1:4
        if tasteResp2(k,t)<0.01
            tasteFlag2(k,t)=1;
        else
            tasteFlag2(k,t)=0;
        end        
    end
end
%% 
% Taste responsivness kruscall wallis
for k = 1:size(Sum,2)   
        Maltose(:,1)    = mean(Sum(k).PSTH.RightRew.Maltose.FRmatrix(:,11:16),2); % take the mean firing rate for each taste
        Maltose(1:size(Maltose,1),2)    = 1; % set the group
        Octacet(:,1)    = mean(Sum(k).PSTH.RightRew.Octacetate.FRmatrix(:,11:16),2);
        Octacet(1:size(Octacet,1),2)    = 2;
        Sucrose(:,1)    = mean(Sum(k).PSTH.LeftRew.Sucrose.FRmatrix(:,11:16),2);
        Sucrose(1:size(Sucrose,1),2)    = 3;
        Quinine(:,1)    = mean(Sum(k).PSTH.LeftRew.Quinine.FRmatrix(:,11:16),2);
        Quinine(1:size(Quinine,1),2)    = 4;
        
        TasteMatrix   = [Maltose;Octacet;Sucrose;Quinine];
        [p(k,1),~,stats]        = kruskalwallis(TasteMatrix(:,1),TasteMatrix(:,2),'off');
        %[c,m,h,gnames] = multcompare(stats);
        if p(k,1)<0.05
              tasteFlag4(k)=1; % the p criteria is 0.05
        else
              tasteFlag4(k)=0;            
        end
        if p(k,1)<0.01
              tasteFlag5(k)=1; % the p criteria is 0.01
        else
              tasteFlag5(k)=0;            
        end
        clear Maltose Sucrose Quinine Octacet
end
%% 
%choose neurons showing both significant response to Kruscall and ranksum
% test.
H=1;
for j=1:214
    if tasteFlag4(1,j)==1 & sum(tasteFlag2(j,:))>0
        TotTasteSel(H,:)=tasteFlag2(j,:);
        Idx(H) = j; % get the id of the taste response
        H = H+1;
    else
    end
end
%% get the taste event for error trials
Sum = tsTasteError(Sum);
%%
% recalculate the psth for each taste response
tasteN = Sum(Idx);
for i = 1:length(tasteN) % reorganize the data
    tastepsth(i).M = tasteN(i).PSTH.RightRew.Maltose;
    tastepsth(i).O = tasteN(i).PSTH.RightRew.Octacetate;
    tastepsth(i).S = tasteN(i).PSTH.LeftRew.Sucrose;
    tastepsth(i).Q = tasteN(i).PSTH.LeftRew.Quinine;
end

taste = {'S', 'M', 'Q','O'};
for i = 1:length(tasteN) % recalculate the psth with 50 ms as bin size, and only take from 0 to 500 ms
    for j = 1:length(taste)
        tastepsth_bin25(i).(taste{j}) = spike2eventRasteandPSTH_NP(tastepsth(i).(taste{j}).Spike, tastepsth(i).(taste{j}).Event, 50,0,500);
    end
    fprintf('Finish processing neuron # %0.f\n',i)
end
%% for error trials for decoding purpose
tasteN = Sum(Idx);
for i = 1:length(tasteN) % reorganize the data
    Errortastepsth(i).M = tasteN(i).event.tsRErr.TasteID.M;
    Errortastepsth(i).O = tasteN(i).event.tsRErr.TasteID.O;
    Errortastepsth(i).S = tasteN(i).event.tsLErr.TasteID.S;
    Errortastepsth(i).Q = tasteN(i).event.tsLErr.TasteID.Q;
end

taste = {'S', 'M', 'Q','O'};
for i = 1:length(tasteN) % recalculate the psth with 50 ms as bin size, and only take from 0 to 500 ms
    for j = 1:length(taste)
        Error_tastepsth_bin25(i).(taste{j}) = spike2eventRasteandPSTH_NP(tasteN(i).timestampN, Errortastepsth(i).(taste{j}), 50,0,500);
    end
    fprintf('Finish processing neuron # %0.f\n',i)
end

%% Normalize the firing rate with auROC
for i = 1:length(tasteN) % recalculate the psth with 25 ms as bin size, and only take from 0 to 500 ms
    for j = 1:length(taste)
        tastepsth_bin25_norm(i).(taste{j}) = psth_auROC(tastepsth_bin25(i).(taste{j}).FRmatrix,0.5, 50);
    end
    fprintf('Finish processing neuron # %0.f\n',i)
end
%%
taste = {'S', 'M', 'Q','O'};
clear X
for i = 1:length(tastepsth_bin25_norm)
    for j = 1:length(taste)
%         X(j,:) = tastepsth_bin25_norm(i).(taste{j})(10:20); % from time 0 to time 0.5s
          X(j,:) = tastepsth_bin25(i).(taste{j}).FR_avg(10:20)
    end
%     D = pdist(X,'correlation');
    D = pdist(X);
    Z = squareform(D);
    DistMatrix(:,:,:,i) = Z;
    DistArray(i,:) = D;
end
Idx_nan = find(isnan((mean(DistArray,2))));
DistMatrix(:,:,:,Idx_nan) = [];
DisMatrix_m = mean(DistMatrix,4);
anova1(DistArray)