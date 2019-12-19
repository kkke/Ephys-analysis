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

%% extract the orofacial data for each session
expe = preLickNoPaw2_reRun(Oro_data,4,1); 
for i = 1:length(expe)
    expe(i).Loro = mean(expe(i).OroFacial.Trial.LSpout,1);
    expe(i).Roro = mean(expe(i).OroFacial.Trial.RSpout,1);   
end
%% plot example trace with licking
i = 5
% figure;
% plot(expe(1).OroFacial.edges,expe(i).Loro)
% hold on
% plot(expe(1).OroFacial.edges,expe(i).Roro)
figure;
plot(Oro_data(1).OroFacial.edges,mean(Oro_data(i).OroFacial.Trial.LSpout))
hold on
plot(Oro_data(1).OroFacial.edges,mean(Oro_data(i).OroFacial.Trial.RSpout))
%%
i=7
K = VideoReader([Oro_data(i).file(1:end-8),'.avi']);
im = read(K,1);
im = im(121:360,211:450,:);
figure;imshow(rgb2gray(im),[0,100])

crop = Oro_data(i).OroFacial.roi(1,:);
hold on

plot(crop(1)-210,crop(2),'o')
plot(crop(1)-210+crop(3),crop(2),'o')
plot(crop(1)-210,crop(2)+crop(4),'o')
plot(crop(1)+crop(3)-210,crop(2)+crop(4),'o')

%%
figure;
plot(expe(7).OroFacial.edges,expe(7).OroFacial.Trial.LSpout(1,:))
xlim([-2,1])
nexData=readNexFile('RVKC235_060217_Final.nex');
time = nexData.contvars{4,1}.timestamps + (0: 1/1000 : (length(nexData.contvars{4,1}.data)-1)*1/1000);
time = time - expe(7).Events.LSpout(1);
idx = find(time>-2 & time<1);
hold on
plot(time(idx),nexData.contvars{4,1}.data(idx))
plot(time(idx),nexData.contvars{3,1}.data(idx))
load('Sum2.mat')
licking = Sum(16).event.raw_decision;
licking_lat = licking{2,8};
figure;
plot(expe(7).OroFacial.edges, expe(7).OroFacial.Trial.LSpout(1,:))
xlim([-2,1])
hold on
scatter(licking_lat-licking_lat(1),ones(size(licking_lat)))
%%
for i = 1:length(expe)
    loro(i,:) = expe(i).Loro;
    roro(i,:) = expe(i).Roro;
end
figure
h1 = boundedline(expe(7).OroFacial.edges, mean(loro,1), std(loro)/sqrt(size(loro,1)),'b')
h2 = boundedline(expe(7).OroFacial.edges, mean(roro,1), std(roro)/sqrt(size(roro,1)),'r')
xlim([-2,1])