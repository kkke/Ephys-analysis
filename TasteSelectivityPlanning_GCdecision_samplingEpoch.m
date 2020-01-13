%% For the planning neurons, let try to extract the information whether they are taste selective
load('Sum2.mat')
clearvars -except Sum
%% Step 1: get the index of the planning neurons
for i = 1:length(Sum)
    resp(i) = Sum(i).auROC.PlanRcorr_vsLcorr.stats.respFlag;
end
idx_p = find(resp ==1);
%% Step2: try to get the psth of maltose  and sucrose octaacetate right trials but aligned to the lateral licks

for j = 1:length(Sum)
    i = j;
    fprintf('Process neuron %4.2f \n',j)
    taste(j).Msampling = spike2eventRasteandPSTH_NP (Sum(i).timestampN, Sum(i).event.tsRCorr.TasteID.M, 500, 0, 500);
    taste(j).Osampling = spike2eventRasteandPSTH_NP (Sum(i).timestampN, Sum(i).event.tsRCorr.TasteID.O, 500, 0, 500);
    [taste(j).MO.p]    = Test_auROC_dISCRIMINATION  (taste(j).Msampling.FRmatrix, taste(j).Osampling.FRmatrix);
    
    % prepare for shuffle to compute statistic
    % put in a single vector the p values from Right and Left
    [taste(j).MO.allp] = [taste(j).Msampling.FRmatrix;taste(j).Osampling.FRmatrix];
    for jj =1:1000
        [taste(j).MO.shuffl.M(:,jj), idx]= datasample(taste(j).MO.allp,size(taste(j).Msampling.FRmatrix,1),'Replace',false);
        idx1 = 1:size(taste(j).MO.allp,1);
        taste(j).MO.shuffl.O(:,jj)       =  taste(j).MO.allp(setdiff(idx1,idx));
        taste(j).MO.shuffl.p (:,jj)       = Test_auROC_dISCRIMINATION(taste(j).MO.shuffl.M(:,jj),taste(j).MO.shuffl.O(:,jj));
    end
    if taste(j).MO.p<0
        [taste(j).MO.stats.pvalue] = length(find(taste(j).MO.shuffl.p<taste(j).MO.p))/1000;
    else
        [taste(j).MO.stats.pvalue] = length(find(taste(j).MO.shuffl.p>taste(j).MO.p))/1000;
    end
    
    if taste(j).MO.stats.pvalue<0.01
        taste(j).MO.stats.respFlag=1;
    else
        taste(j).MO.stats.respFlag=0;
    end
end
%% Step3: try to get the psth of sucrose  and  quinine left trials but aligned to the lateral licks

for j = 1:length(Sum)
    fprintf('Process neuron %4.2f \n',j)
    i = j;
    taste(j).Splanning = spike2eventRasteandPSTH_NP (Sum(i).timestampN, Sum(i).event.tsLCorr.TasteID.S, 500, 0, 500);
    taste(j).Qplanning = spike2eventRasteandPSTH_NP (Sum(i).timestampN,Sum(i).event.tsLCorr.TasteID.Q, 500, 0, 500);
    [taste(j).SQ.p]    = Test_auROC_dISCRIMINATION  (taste(j).Splanning.FRmatrix, taste(j).Qplanning.FRmatrix);
    
    % prepare for shuffle to compute statistic
    % put in a single vector the p values from Right and Left
    [taste(j).SQ.allp] = [taste(j).Splanning.FRmatrix;taste(j).Qplanning.FRmatrix];
    for jj =1:1000
        [taste(j).SQ.shuffl.S(:,jj), idx]= datasample(taste(j).SQ.allp,size(taste(j).Splanning.FRmatrix,1),'Replace',false);
        idx1 = 1:size(taste(j).SQ.allp,1);
        taste(j).SQ.shuffl.Q(:,jj)       =  taste(j).SQ.allp(setdiff(idx1,idx));
        taste(j).SQ.shuffl.p (:,jj)       = Test_auROC_dISCRIMINATION(taste(j).SQ.shuffl.S(:,jj),taste(j).SQ.shuffl.Q(:,jj));
    end
    if taste(j).SQ.p<0
        [taste(j).SQ.stats.pvalue] = length(find(taste(j).SQ.shuffl.p<taste(j).SQ.p))/1000;
    else
        [taste(j).SQ.stats.pvalue] = length(find(taste(j).SQ.shuffl.p>taste(j).SQ.p))/1000;
    end
    
    if taste(j).SQ.stats.pvalue<0.01
        taste(j).SQ.stats.respFlag=1;
    else
        taste(j).SQ.stats.respFlag=0;
    end
end
%% try to summarize to see which neurons show taste selectivity
for i = 1:length(taste)
    fprintf('Process neuron %4.2f \n',i)
    if taste(i).MO.stats.pvalue ==0 
       shu = taste(i).MO.shuffl.p;
       ok_p = length(find(taste(i).MO.p >= shu))/1000;
       oko_p= length(find(taste(i).MO.p <= shu))/1000;
       if ok_p >=0.01 & oko_p >=0.01 % if this exist, it means shuffling preference and real preference has lots in common
           fprintf('Neuron %4.2f is not good Left\n',i)
           taste(i).MO.stats.respFlag = 0;
       end
    end
       
    if taste(i).SQ.stats.pvalue ==0 
       shu2 = taste(i).SQ.shuffl.p;
       ok_p = length(find(taste(i).SQ.p >= shu2))/1000;
       oko_p= length(find(taste(i).SQ.p <= shu2))/1000;
       if ok_p >=0.01 & oko_p >=0.01 % if this exist, it means shuffling preference and real preference has lots in common
           fprintf('Neuron %4.2f is not good Right\n',i)
           taste(i).SQ.stats.respFlag = 0;
       end
    end
end
%%
%% Step4: try to get the psth of left trials but aligned to the central licks
for j = 1:length(Sum)
    fprintf('Process neuron %4.2f \n',j)
    i = j;
    event_L = sort([Sum(i).event.tsLCorr.TasteID.S; Sum(i).event.tsLCorr.TasteID.Q]);
    taste(j).SQplanning = spike2eventRasteandPSTH_NP (Sum(i).timestampN, event_L, 500, 0, 500);
    event_R = sort([Sum(i).event.tsRCorr.TasteID.M; Sum(i).event.tsRCorr.TasteID.O]);
    taste(j).MOplanning = spike2eventRasteandPSTH_NP (Sum(i).timestampN,event_R, 500, 0, 500);
    [taste(j).LR.p]    = Test_auROC_dISCRIMINATION  (taste(j).SQplanning.FRmatrix, taste(j).MOplanning.FRmatrix);
    
    % prepare for shuffle to compute statistic
    % put in a single vector the p values from Right and Left
    [taste(j).LR.allp] = [taste(j).SQplanning.FRmatrix;taste(j).MOplanning.FRmatrix];
    for jj =1:1000
        [taste(j).LR.shuffl.L(:,jj), idx]= datasample(taste(j).LR.allp,size(taste(j).SQplanning.FRmatrix,1),'Replace',false);
        idx1 = 1:size(taste(j).LR.allp,1);
        taste(j).LR.shuffl.R(:,jj)       =  taste(j).LR.allp(setdiff(idx1,idx));
        taste(j).LR.shuffl.p(:,jj)       = Test_auROC_dISCRIMINATION(taste(j).LR.shuffl.L(:,jj),taste(j).LR.shuffl.R(:,jj));
    end
    if taste(j).LR.p<0
        [taste(j).LR.stats.pvalue] = length(find(taste(j).LR.shuffl.p<taste(j).LR.p))/1000;
    else
        [taste(j).LR.stats.pvalue] = length(find(taste(j).LR.shuffl.p>taste(j).LR.p))/1000;
    end
    
    if taste(j).LR.stats.pvalue<0.01
        taste(j).LR.stats.respFlag=1;
    else
        taste(j).LR.stats.respFlag=0;
    end
end
%%
%% try to summarize to see which neurons show taste selectivity
for i = 1:length(taste)
    fprintf('Process neuron %4.2f \n',i)
    if taste(i).LR.stats.pvalue ==0 
       shu = taste(i).LR.shuffl.p;
       ok_p = length(find(taste(i).LR.p >= shu))/1000;
       oko_p= length(find(taste(i).LR.p <= shu))/1000;
       if ok_p >=0.01 & oko_p >=0.01 % if this exist, it means shuffling preference and real preference has lots in common
           fprintf('Neuron %4.2f is not good Left\n',i)
           taste(i).LR.stats.respFlag = 0;
       end     
    else
    end
end
save('Planning_Taste_Selectivity_allNeuron_samplingEpoch.mat','taste')
%%
clear bothnon
x=1;
for i = 1:length(taste)
       if taste(i).MO.stats.respFlag || taste(i).SQ.stats.respFlag
           fprintf('Neuron %4.2f is taste selective \n',i)
           bothnon(x) = i;
           x = x+1;
       end
end
clear bothnon_plan
x=1;
for i = 1:length(taste)
       if taste(i).LR.stats.respFlag
           fprintf('Neuron %4.2f is taste selective \n',i)
           bothnon_plan(x) = i;
           x = x+1;
       end
end

%% plot the left/right preference vs taste preference
figure;
for i = 1:length(idx_p) % plot the planning and taste selectivity during the sampling epoch.
    taste_pall(i) = max([abs(taste(idx_p(i)).MO.p),abs(taste(idx_p(i)).SQ.p)]);
    lr_pall(i)    = taste(idx_p(i)).LR.p;
end
plot(abs(taste_pall),abs(lr_pall),'bo')

bothnon1 = intersect(bothnon,idx_p); % taste selectivity(sample)+ planning (delay)
bothnon2 = intersect(bothnon_plan,idx_p); % planning (sample) + planning (delay)
for i = 1:length(bothnon1)
    taste_p(i) = max([abs(taste(bothnon1(i)).MO.p),abs(taste(bothnon1(i)).SQ.p)]);
    lr_p(i)    = taste(bothnon1(i)).LR.p;
end

% for i = 1:length(bothnon2) % for highlighting the still planning
%     taste_pp(i) = max([abs(taste(bothnon2(i)).MO.p),abs(taste(bothnon2(i)).SQ.p)]);
%     lr_pp(i)    = taste(bothnon2(i)).LR.p;
% end
% hold on 
% plot(abs(taste_pp),abs(lr_pp),'go')

hold on 
plot(abs(taste_p),abs(lr_p),'ro')
xlim([0,1])
ylim([0,1])
hold on
plot([0,1],[0,1])
axis square
xlabel('Max taste selectivity')
ylabel('Planning selectivity')
%%
save('PlanningTaste.mat','taste')
%% try to show neurons
% close all

i = bothnon(20); % good example: bothnon(29) for high/high; i = 33,49, bothnon(10,34),for High/Low, bothnon(2),bothnon(3), bothnon(23) for low/high; bothnon(7,17,32) for low/low
M = spike2eventRasteandPSTH_NP(taste(i).Mplanning.Spike,taste(i).Mplanning.Event,100,-2000,1000);
O = spike2eventRasteandPSTH_NP(taste(i).Oplanning.Spike,taste(i).Oplanning.Event,100,-2000,1000);
S = spike2eventRasteandPSTH_NP(taste(i).Splanning.Spike,taste(i).Splanning.Event,100,-2000,1000);
Q = spike2eventRasteandPSTH_NP(taste(i).Qplanning.Spike,taste(i).Qplanning.Event,100,-2000,1000);

% figure;
% figure; plot(a.timepoint,a.FR_avg)
% hold on
% plot(b.timepoint,b.FR_avg)
% plot(c.timepoint,c.FR_avg)
% plot(d.timepoint,d.FR_avg)
% legend({'M','O','S','Q'})
% xlim([-2,1])

    %%%%%%%%%%%%%%-----------Figure part- -----------%%%%%%%%%%%%%%%%
    
    colors=distinguishable_colors(5); %generate some vectors contaning col values
    
    figure;
    %%raster
    subplot (2,1,1);
    tot_tr=[size(S.spikeraster,2), size(Q.spikeraster,2), size(M.spikeraster,2), size(O.spikeraster,2)];
    
    for t = 1:sum(tot_tr)
        if (t <= (tot_tr(1)))
            for sp = 1:size(S.spikeraster(t).times,2)
                %                 Lin = line([S.spikeraster(t).times(sp) S.spikeraster(t).times(sp)],...
                %                     [t-1 t]);
                Lin = plot([S.spikeraster(t).times(sp),S.spikeraster(t).times(sp)],[t-1 t])
                hold on
                set(Lin,'Color',colors(1,:),'LineWidth',1);
            end
        elseif (t > (tot_tr(1)) &&  t <= (tot_tr(1)+tot_tr(2)))
            for sp = 1:size(Q.spikeraster(t-size(S.spikeraster,2)).times,2)
                Lin1 = plot([Q.spikeraster((t-size(S.spikeraster,2))).times(sp) Q.spikeraster((t-size(S.spikeraster,2))).times(sp)],...
                    [t-1 t]);
                set(Lin1,'Color',colors(2,:),'LineWidth',1);
            end
        elseif (t > (tot_tr(1)+tot_tr(2))&& t <= (tot_tr(1)+tot_tr(2)+tot_tr(3)))
            for sp = 1:size(M.spikeraster(t-(size(S.spikeraster,2)+size(Q.spikeraster,2))).times,2)
                Lin1 = plot([M.spikeraster(t-(size(S.spikeraster,2)+size(Q.spikeraster,2))).times(sp) M.spikeraster(t-(size(S.spikeraster,2)+size(Q.spikeraster,2))).times(sp)],...
                    [t-1 t]);
                set(Lin1,'Color',colors(3,:),'LineWidth',1);
            end
        elseif (t > (tot_tr(1)+tot_tr(2)+tot_tr(3))&& t <= (tot_tr(1)+tot_tr(2)+tot_tr(3)+tot_tr(4)))
            for sp = 1:size(O.spikeraster(t-(size(S.spikeraster,2)+size(Q.spikeraster,2)+size(M.spikeraster,2))).times,2)
                Lin1 = plot([O.spikeraster(t-(size(S.spikeraster,2)+size(Q.spikeraster,2)+size(M.spikeraster,2))).times(sp) O.spikeraster(t-(size(S.spikeraster,2)+size(Q.spikeraster,2)+size(M.spikeraster,2))).times(sp)],...
                    [t-1 t]);
                set(Lin1,'Color',colors(4,:),'LineWidth',1);
            end
        end
    end
    hold on;
    %line([0 0],[0 sum(tot_tr)],'LineWidth',0.001,'Color','k');
    box('off') ;
    axis('off');
    ylabel(gca,'Trials');
    set(get(gca,'YLabel'),'Visible','on','fontsize',12);
%     title([forFig.mouse{n,1} ' _ ' forFig.date{n,1} ' _ ' forFig.neuron{n,1} ' _ bin=' num2str(bin) 'ms']);
    %
    % psth
    tsize = 3;
    tsig  = 2;
%     S.FR_avg = gaussmooth(S.FR_avg,tsize,tsig);
%     Q.FR_avg = gaussmooth(Q.FR_avg,tsize,tsig);
%     M.FR_avg = gaussmooth(M.FR_avg,tsize,tsig);
%     O.FR_avg = gaussmooth(O.FR_avg,tsize,tsig);
    subplot (2,1,2);
    h = plot(S.timepoint,S.FR_avg);hold on;
    set(h,'Color',colors(1,:));
    h1= plot(S.timepoint,Q.FR_avg);hold on;
    set(h1,'Color',colors(2,:));
    h2= plot(S.timepoint,M.FR_avg);hold on;
    set(h2,'Color',colors(3,:));
    h3= plot(S.timepoint,O.FR_avg);hold on;
    set(h3,'Color',colors(4,:));
    %h3= plot([0 0],[0 round(max([R.FR_avg L.FR_avg]))],'--k');
    
    % set x axis
%     xlim([time(1) time(2)]);
%     xlabel('Time(s)');
%     xticks([time(1) time(2)]);
%     xticklabels({num2str(time(1)),num2str(time(2))});
    
    % set y axis
    % set y axis
    ylim([0 10]);
%     yticks([0 14]);
%     yticklabels({'0','14'});
    ylabel('Firing rate');
    
    %title('goal-selective neurons')
    set(gca,'fontsize',12);
    box('off') ;
    legend('S','Q','M','O');
    
    % prepare for annotation
    xa = [0 0]; % X-Coordinates in data space
    ya = [0 sum(tot_tr)]; % Y-Coordinates in data space
    [xaf,yaf] = ds2nfu(xa,ya);
    yaf(2)=0.9;
    annotation('line',xaf,yaf);
    clearvars -except Sum taste taste_p bothnon idx_p lr_p lr_pall taste_pall;