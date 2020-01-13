%%
load('LickingParameter_Chemogenetics.mat')
%% summarize the sampling duration and reaction time
control.saline = licking_param_chemogenetic(control.saline);
control.cno    = licking_param_chemogenetic(control.cno);
hM4Di.saline = licking_param_chemogenetic(hM4Di.saline);
hM4Di.cno    = licking_param_chemogenetic(hM4Di.cno);
