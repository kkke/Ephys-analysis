%function fileForFig_rastPSTH
% example figure for FYAP grants
% licking behavior 5- 6 licks before a taste is delivered


% close all

n   =18;

tsN = forFig.timestampN{n,1};
%time around the event you want to use for the plot
time = [-2 1];
%binsize
bin=100;

   
% for correct trials
%     tsEv_1 = forFig.event{n,1}.tsRCorr.FLickRSpou;
%     tsEv_2 = forFig.event{n,1}.tsLCorr.FLickLSpou;
% for error  trials:
    tsEv_1 = forFig.event{n,1}.tsRErr.FLickRSpou;
    tsEv_2 = forFig.event{n,1}.tsLErr.FLickRSpou;

    %extracting the PSTH
    R = spike2eventRasteandPSTH_NP (tsN,tsEv_1, bin, time(1)*1000, time(2)*1000);
    L = spike2eventRasteandPSTH_NP (tsN,tsEv_2, bin, time(1)*1000, time(2)*1000);
    
    %%%%%%%%%%%%%%-----------Figure part- -----------%%%%%%%%%%%%%%%%
    
    colors=distinguishable_colors(5); %generate some vectors contaning col values
    
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
    title([forFig.mouse{n,1} ' _ ' forFig.date{n,1} ' _ ' forFig.neuron{n,1} ' _ bin=' num2str(bin) 'ms']);
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
ylim([0,32])
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
%     clearvars -except forFig Sum;
    %%
    %         for t = 1:size(data{1,s}.psth{j,1}.spikeraster,2)
    %             for sp = 1:size(data{1,s}.psth{j,1}.spikeraster(t).times,2)
    %                 line([data{1,s}.psth{j,1}.spikeraster(t).times(sp) data{1,s}.psth{j,1}.spikeraster(t).times(sp)],...
    %                     [t-1 t],'Color','k');
    %             end
    %         end
    
  