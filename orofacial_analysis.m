filelist = dir('*Oro.mat');
for i = 1:length(filelist)
    load(filelist(i).name)
    Oro_data(i).file = filelist(i).name;
    Oro_data(i).OroFacial = expe.OroFacial;
    Oro_data(i).Events    = expe.Events;
    Oro_data(i).Times     = expe.Times;
    Oro_data(i).auROC     = expe.auROC;   
end

for i = 1:length(filelist)
    Rspout(i,:) = Oro_data(i).auROC.RSpout';
    Lspout(i,:) = Oro_data(i).auROC.LSpout';
end

figure;
h1 = boundedline(Oro_data(1).OroFacial.edges, mean(Rspout,1), std(Rspout)/sqrt(size(Rspout,1)),'r')
h2 = boundedline(Oro_data(1).OroFacial.edges, mean(Lspout,1), std(Rspout)/sqrt(size(Lspout,1)),'b')

%%
expe = preLickNoPaw2_reRun(Oro_data,4,1) 