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
%%
