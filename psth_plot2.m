%function fileForFig_rastPSTH
% example figure for FYAP grants
% licking behavior 5- 6 licks before a taste is delivered


% close all
function psth_plot2(forFig, n,c)
% n   =18;
eventNo = 2;
tsN = forFig(n).timestampN;

%time around the event you want to use for the plot
time = [-2 1];
% time = [-.5 0.5];

%binsize
 bin=100; % for planning/goal
%  bin=25;% for taste

% current neurons for representative figure are:
% taste: 13, 19, 27,34,39,43, 70, 77, 96
% goal: neuron#1 is 96, neuron# is 33
% 18,91,102 is the first example for planning neuron; to plot also the orofacial movement of
% this session load('RVKC235_060217_Oro.mat') and run:
%boundedline(expe.OroFacial.edges, mean(expe.OroFacial.Trial.RSpout(:,:),1), std(expe.OroFacial.Trial.RSpout(:,:),1)./sqrt(size(expe.OroFacial.Trial.RSpout,1)),'r');hold on;
%boundedline(expe.OroFacial.edges, mean(expe.OroFacial.Trial.LSpout(:,:),1), std(expe.OroFacial.Trial.LSpout(:,:),1)./sqrt(size(expe.OroFacial.Trial.LSpout,1)),'b');hold on;
% 112 is the second example for planning neuron; to plot also the orofacial movement of
% this session load('RVKC235_053117_Oro.mat') and run:
% 79,96,115,127,134 is the one potential example for goal neuron;
% for representative neuron tracking the planning: 18, 113, 134,155,158,202,208,212,214 ;
% for representative neurons not trakcing the planning: 48,49, 106, 139,148
switch eventNo
case 2
    %event
    if c == 1
        tsEv_1 = forFig(n).event.tsRCorr.FLickRSpou;
        tsEv_2 = forFig(n).event.tsLCorr.FLickLSpou;
%         tsEv_2 = forFig(n).event.tsRErr.FLickRSpou;
        
    else
        % for error  trials:
        tsEv_1 = forFig(n).event.tsRErr.FLickRSpou;
        tsEv_2 = forFig(n).event.tsLErr.FLickRSpou;
%         tsEv_2 = forFig(n).event.tsLCorr.FLickLSpou;
        
    end

    %extracting the PSTH
    R = spike2eventRasteandPSTH_NP (tsN,tsEv_1, bin, time(1)*1000, time(2)*1000);
    L = spike2eventRasteandPSTH_NP (tsN,tsEv_2, bin, time(1)*1000, time(2)*1000);
    
    %%%%%%%%%%%%%%-----------Figure part- -----------%%%%%%%%%%%%%%%%
    
    colors=distinguishable_colors(5); %generate some vectors contaning col values
    
    figure
%     figure(n);
    %%raster
    subplot (2,1,1);
    tot_tr=[size(R.spikeraster,2) size(L.spikeraster,2)];
    
    for t = 1:sum(tot_tr)
        if (t < (tot_tr(1))) || (t == (tot_tr(1)))
            for sp = 1:size(R.spikeraster(t).times,2)
                Lin = line([R.spikeraster(t).times(sp) R.spikeraster(t).times(sp)],...
                    [t-1 t]);
                %Lin = plot([R.spikeraster(t).times(sp),R.spikeraster(t).times(sp)],[t-1 t]);
                hold on    
                set(Lin,'Color',colors(1,:),'LineWidth',1);
                set(Lin,'Color',colors(1,:),'LineWidth',1);
            end
        else
            for sp = 1:size(L.spikeraster(t-size(R.spikeraster,2)).times,2)
                Lin1 = line([L.spikeraster((t-size(R.spikeraster,2))).times(sp) L.spikeraster((t-size(R.spikeraster,2))).times(sp)],...
                    [t-1 t]);
                %Lin1 = plot([L.spikeraster(t).times(sp),L.spikeraster(t).times(sp)],[t-1 t]);

                set(Lin1,'Color',colors(3,:),'LineWidth',1);
            end
        end
    end
    hold on;
    %line([0 0],[0 sum(tot_tr)],'LineWidth',0.001,'Color','k');
    box('off') ;
    axis('off');
    ylabel(gca,'Trials');
    set(get(gca,'YLabel'),'Visible','on','fontsize',12);
    title([forFig(n).mouse ' _ ' forFig(n).date ' _ ' forFig(n).neuron ' _ bin=' num2str(bin) 'ms']);
    %
    % psth
    subplot (2,1,2);
    h = plot(R.timepoint,R.FR_avg);hold on;
    set(h,'Color',colors(1,:));
    h1= plot(L.timepoint,L.FR_avg);hold on;
    set(h1,'Color',colors(3,:));
    %h3= plot([0 0],[0 round(max([R.FR_avg L.FR_avg]))],'--k');
    
    % set x axis
    xlim([time(1) time(2)]);
    xlabel('Time(s)');
    xticks([time(1) time(2)]);
    xticklabels({num2str(time(1)),num2str(time(2))});
    
    % set y axis
%     ylim([0 round(max([R.FR_avg L.FR_avg]))]);
%     yticks([0 round(max([R.FR_avg L.FR_avg]))]);
%     yticklabels({'0',num2str(round(max([R.FR_avg L.FR_avg])))});
% ylim([0,32])
    ylabel('Firing rate');
    
    %title('goal-selective neurons')
    set(gca,'fontsize',12);
    box('off') ;
    legend('Right-Corr','Left-Corr');
    
    % prepare for annotation
    xa = [0 0]; % X-Coordinates in data space
    ya = [0 sum(tot_tr)]; % Y-Coordinates in data space
    [xaf,yaf] = ds2nfu(xa,ya);
    yaf(2)=0.9;
    annotation('line',xaf,yaf);
    clearvars -except forFig Sum;
    %%
    %         for t = 1:size(data{1,s}.psth{j,1}.spikeraster,2)
    %             for sp = 1:size(data{1,s}.psth{j,1}.spikeraster(t).times,2)
    %                 line([data{1,s}.psth{j,1}.spikeraster(t).times(sp) data{1,s}.psth{j,1}.spikeraster(t).times(sp)],...
    %                     [t-1 t],'Color','k');
    %             end
    %         end
    
    case 4
    tsEv_1 = forFig.event{n, 1}.tsLCorr.TasteID.S;
    tsEv_2 = forFig.event{n, 1}.tsLCorr.TasteID.Q;
    tsEv_3 = forFig.event{n, 1}.tsRCorr.TasteID.M;
    tsEv_4 = forFig.event{n, 1}.tsRCorr.TasteID.O;
    %extracting the PSTH
%     S = spike2eventRasteandPSTH_NP (SessionData.PlexonTsNeuron.Unit1.ts, SessionData.PlexonTsTaste.Sucrose, 100, -1000, 1000);
%     Q = spike2eventRasteandPSTH_NP (SessionData.PlexonTsNeuron.Unit1.ts, SessionData.PlexonTsTaste.Quinine, 100, -1000, 1000)
%     M = spike2eventRasteandPSTH_NP (SessionData.PlexonTsNeuron.Unit1.ts, SessionData.PlexonTsTaste.NaCl, 100, -1000, 1000)
%     O = spike2eventRasteandPSTH_NP (SessionData.PlexonTsNeuron.Unit1.ts, SessionData.PlexonTsTaste.NaCl, 100, -1000, 1000)
    S = spike2eventRasteandPSTH_NP (forFig.timestampN{n,1}, tsEv_1, bin, time(1)*1000, time(2)*1000);
    Q = spike2eventRasteandPSTH_NP (forFig.timestampN{n,1}, tsEv_2, bin, time(1)*1000, time(2)*1000);
    M = spike2eventRasteandPSTH_NP (forFig.timestampN{n,1}, tsEv_3, bin, time(1)*1000, time(2)*1000);
    O = spike2eventRasteandPSTH_NP (forFig.timestampN{n,1}, tsEv_4, bin, time(1)*1000, time(2)*1000);
    %%%%%%%%%%%%%%-----------Figure part- -----------%%%%%%%%%%%%%%%%
    
    colors=distinguishable_colors(5); %generate some vectors contaning col values
    
    figure(n);
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
    title([forFig.mouse{n,1} ' _ ' forFig.date{n,1} ' _ ' forFig.neuron{n,1} ' _ bin=' num2str(bin) 'ms']);
    %
    % psth
    tsize = 3;
    tsig  = 2;
    S.FR_avg = gaussmooth(S.FR_avg,tsize,tsig);
    Q.FR_avg = gaussmooth(Q.FR_avg,tsize,tsig);
    M.FR_avg = gaussmooth(M.FR_avg,tsize,tsig);
    O.FR_avg = gaussmooth(O.FR_avg,tsize,tsig);
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
    xlim([time(1) time(2)]);
    xlabel('Time(s)');
    xticks([time(1) time(2)]);
    xticklabels({num2str(time(1)),num2str(time(2))});
    
    % set y axis
    ylim([0 round(max([S.FR_avg M.FR_avg O.FR_avg Q.FR_avg]))+5]);
    yticks([0 round(max([S.FR_avg M.FR_avg O.FR_avg Q.FR_avg]))+5]);
    yticklabels({'0',num2str(round(max([S.FR_avg M.FR_avg O.FR_avg Q.FR_avg])))});
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
    clearvars -except forFig TasteResID;
end
        