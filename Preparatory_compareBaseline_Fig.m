clear baseline L_firing_r R_firing_r
for i = 1:length(idx_p)
%     baseline(i)  = 1/2*mean(taste(idx_p_R(i)).Lplanning.FR_avg(1:2) + taste(idx_p_R(i)).Rplanning.FR_avg(1:2));
    L_firing_r(i)  = Sum(idx_p(i)).auROC.PlanRcorr_vsLcorr.Left.FR_avg-mean(taste(idx_p(i)).Lplanning.FR_avg(1:2));
    R_firing_r(i)  = Sum(idx_p(i)).auROC.PlanRcorr_vsLcorr.Right.FR_avg-mean(taste(idx_p(i)).Rplanning.FR_avg(1:2));
end
indx = find(abs(L_firing_r)>abs(R_firing_r));
indx_r = setdiff(1:length(idx_p),indx);
for i = 1:length(indx)
    p_l(i) = Sum(idx_p(indx(i))).auROC.PlanRcorr_vsLcorr.p;

end

for i = 1:length(indx_r)
    p_r(i) = Sum(idx_p(indx_r(i))).auROC.PlanRcorr_vsLcorr.p;

end
non_idx = setdiff(1:length(Sum),idx_p);
for i = 1:length(non_idx)
    p_n(i) = Sum(non_idx(i)).auROC.PlanRcorr_vsLcorr.p;

end


figure
h = histogram(p_l,30);hold on;
h1= histogram(p_r,30);
h3 = histogram(p_n,30);hold on;
bin      = 0.06;
%histogram color settings
colors=distinguishable_colors(5);

set(h,'FaceColor',colors(1,:),'EdgeColor','w','BinWidth',bin);
set(h1,'FaceColor',colors(3,:),'EdgeColor','w','BinWidth',bin);
set(h3,'FaceColor',colors(2,:),'EdgeColor','w','BinWidth',bin);

xlim([-1 1]);
xticks([-1 -0.5 0 0.5 1]);
yticks([0 30]);
ylabel('Neurons');
yticklabels({'0',num2str(30)});
ylim([0,30])
