%function fileForFig_rastPSTH
% example figure for FYAP grants
% licking behavior 5- 6 licks before a taste is delivered


% close all
% n = 1, 7, 13

n   =13;

%time around the event you want to use for the plot
time = [-1 1];
%binsize
bin=50;

figure
   
% for correct trials
    R   = spike2eventRasteandPSTH_NP(Sum(n).timestampN,outcome(n).dscorr,bin, time(1)*1000, time(2)*1000);
    % error
    L   = spike2eventRasteandPSTH_NP(Sum(n).timestampN,outcome(n).error,bin, time(1)*1000, time(2)*1000);
% for licking
%     R   = spike2eventRasteandPSTH_NP(outcome(i).corLicks,outcome(i).dscorr,bin, time(1)*1000, time(2)*1000);
%     L = spike2eventRasteandPSTH_NP(outcome(i).errLicks,outcome(i).error,bin, time(1)*1000, time(2)*1000);
    
    %extracting the PSTH

    
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
    xlim([-1,1])
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
%     xticks([time(1) time(2)]);
%     xticklabels({num2str(time(1)),num2str(time(2))});
    
    % set y axis
%     ylim([0 round(max([R.FR_avg L.FR_avg]))]);
%     yticks([0 round(max([R.FR_avg L.FR_avg]))]);
%     yticklabels({'0',num2str(round(max([R.FR_avg L.FR_avg])))});
    ylim([0,30])
    ylabel('Firing rate');
    
    %title('goal-selective neurons')
    set(gca,'fontsize',12);
    box('off') ;
    legend('Correct','Error');
    
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
    
  