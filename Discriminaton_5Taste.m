%% behaviral analysis for 2p rig
%% Step 1: load data from Intan
% add directory
clear
addpath('C:\Users\ke-roberto\Documents\MATLAB\Code_Ke\')
%% load data
file = dir('*.rhd');
dataRaw = read_Intan(file.name);
%% extract the events
thr = 2;
[data.centSp.dig,~] = Timing_onset_offset(dataRaw.analog(1,:), dataRaw.ts, thr,30,0); % get the central licks
[data.LeftSp.dig,~] = Timing_onset_offset(dataRaw.analog(2,:), dataRaw.ts, thr,30,0);
[data.RightSp.dig,~] = Timing_onset_offset(dataRaw.analog(3,:), dataRaw.ts, thr,30,0);
[data.Up.dig,data.Up.dig_offset]                        = Timing_onset_offset(dataRaw.event(10,:), dataRaw.ts, 0.5,30,0);
[data.Forward.dig,data.Forward.dig_offset]              = Timing_onset_offset(dataRaw.event(9,:), dataRaw.ts, 0.5,30,0);
[data.Down.dig,data.Down.dig_offset]                    = Timing_onset_offset(dataRaw.event(13,:), dataRaw.ts, 0.5,30,0);
%%
%% events
A = cd;
data.mouseID  = A(45:end);
data.date     = A(34:43);
clear A;
%%
[R_1, R_1_2]            = Timing_onset_offset(dataRaw.event(8,:), dataRaw.ts, 0.5,30,0);
[L_1, L_1_2]            = Timing_onset_offset(dataRaw.event(7,:), dataRaw.ts, 0.5,30,0);
[Taste_5,Taste_5_2]     = Timing_onset_offset(dataRaw.event(5,:), dataRaw.ts, 0.5,30,0);
[Taste_4,Taste_4_2]     = Timing_onset_offset(dataRaw.event(4,:), dataRaw.ts, 0.5,30,0);
[Taste_3,Taste_3_2]     = Timing_onset_offset(dataRaw.event(3,:), dataRaw.ts, 0.5,30,0);
[Taste_2,Taste_2_2]     = Timing_onset_offset(dataRaw.event(2,:), dataRaw.ts, 0.5,30,0);
[Taste_1,Taste_1_2]     = Timing_onset_offset(dataRaw.event(1,:), dataRaw.ts, 0.5,30,0);
Taste_1 = sort([Taste_1,Taste_1_2]);
Taste_2 = sort([Taste_2,Taste_2_2]);
Taste_3 = sort([Taste_3,Taste_3_2]);
Taste_4 = sort([Taste_4,Taste_4_2]);
Taste_5 = sort([Taste_5,Taste_5_2]);
L_1     = sort([L_1,L_1_2]);
R_1     = sort([R_1,R_1_2]);

%% reorganize the data
if isempty(Taste_1)
    Taste_1 =[];
else
    Taste_1(2,:) = 0;
end

if isempty(Taste_2)
    Taste_2 =[];
else
    Taste_2(2,:) = 1;
end

if isempty(Taste_3)
    Taste_3 =[];
else
    Taste_3(2,:) = 2;
end

if isempty(Taste_4)
    Taste_4 =[];
else
    Taste_4(2,:) = 3;
end

if isempty(Taste_5)
    Taste_5 =[];
else
    Taste_5(2,:) = 4;
end

if isempty(L_1)
    L_1 =[];
else
    L_1(2,:) = 6;
end

if isempty(R_1)
    R_1 =[];
else
    R_1(2,:) = 7;
end

Taste= [Taste_1, Taste_2, Taste_3, Taste_4, Taste_5,L_1,R_1];
[timestamps1,Idx] = sort(Taste(1,:));
for i = 1:length(Idx)
    data1(i) = Taste(2,Idx(i));
end
timestamps1 = timestamps1';
%%
%%
%data1(1:2)                 =[]; % remove the first two values cause they are useless.
%timestamps1(1:2)           =[]; % remove the first two values cause they are useless.
Left_Correct               = timestamps1(data1==6);
data.Left_Correct.dig      = Left_Correct(1:2:length(Left_Correct));
Right_Correct               = timestamps1(data1==7);
data.Right_Correct.dig      = Right_Correct(1:2:length(Right_Correct));
a=unique(data1);                % get the id of stimuli
tasteID=a(find(a<6));
message1=['You train the animal with ',num2str(length(tasteID)),' tastant, Please specify the tastant for each line.'];
message2=[num2str((tasteID(:)+1)')];
uiwait(msgbox({message1, message2}));
x=input('Specify the tastant\n');
y=input('Specify the left-right of each tastant,1 is left and 2 is right\n')
switch length(tasteID)
        case 2
        tastant_1                                                =timestamps1(data1==a(1));
        data.tastant_1.dig                                       =tastant_1(1:2:length(tastant_1));
%         [data.tastant_1.psth_raster] = spike2eventRasteandPSTH_NP (data.centSp.dig, data.tastant_1.dig, 100, -5000, 5000);
        tastant_2                                                =timestamps1(data1==a(2));
        data.tastant_2.dig                                       =tastant_2(1:2:length(tastant_2 ));
%         [data.tastant_2.psth_raster] = spike2eventRasteandPSTH_NP (data.centSp.dig, data.tastant_2.dig, 100, -5000, 5000);
        data.tastant_1.id=x{1};
        data.tastant_1.line=tasteID(1)+1;
        data.tastant_1.LR=y(1);
        data.tastant_2.id=x{2};
        data.tastant_2.line=tasteID(2)+1;
        data.tastant_2.LR=y(2);
        if find(y==1)==1
            data.leftID='tastant_1';
            data.rightID='tastant_2';
        else
            data.leftID='tastant_2';
            data.rightID='tastant_1';
        end
%         for i=1:length(data.tastant_1.psth_raster.spikeraster)
%             h1=scatter(data.tastant_1.psth_raster.spikeraster(i).times, i*ones(size(data.tastant_1.psth_raster.spikeraster(i).times)),6,'cv','filled');
%             hold on
%         end
%         for j=1:length(data.tastant_2.psth_raster.spikeraster)
%             i=i+1;
%             h2=scatter(data.tastant_2.psth_raster.spikeraster(j).times, i*ones(size(data.tastant_2.psth_raster.spikeraster(j).times)),6,'mv','filled');
%             hold on
%         end
%         legend([h1,h2],{x{1},x{2}})
        
%     case 3
%         tastant_1                                                =timestamps1(data1==a(1));
%         data.tastant_1.dig                                       =tastant_1(1:2:length(tastant_1));
% %         [data.tastant_1.psth_raster] = spike2eventRasteandPSTH_NP (data.centSp.dig, data.tastant_1.dig, 100, -5000, 5000);
%         tastant_2                                                =timestamps1(data1==a(2));
%         data.tastant_2.dig                                       =tastant_2(1:2:length(tastant_2 ));
% %         [data.tastant_2.psth_raster] = spike2eventRasteandPSTH_NP (data.centSp.dig, data.tastant_2.dig, 100, -5000, 5000);
%         tastant_3                                                =timestamps1(data1==a(3));
%         data.tastant_3.dig                                       =tastant_3(1:2:length(tastant_3 ));
% %         [data.tastant_3.psth_raster] = spike2eventRasteandPSTH_NP (data.centSp.dig, data.tastant_3.dig, 100, -5000, 5000);
%         data.tastant_1.id=x{1};
%         data.tastant_1.LR=y(1);
%         data.tastant_2.id=x{2};
%         data.tastant_2.LR=y(2);
%         data.tastant_3.id=x{3};
%         data.tastant_3.LR=y(3);        
%         for i=1:length(data.tastant_1.psth_raster.spikeraster)
%             h1=scatter(data.tastant_1.psth_raster.spikeraster(i).times, i*ones(size(data.tastant_1.psth_raster.spikeraster(i).times)),6,'cv','filled');
%             hold on
%         end
%         for j=1:length(data.tastant_2.psth_raster.spikeraster)
%             i=i+1;
%             h2=scatter(data.tastant_2.psth_raster.spikeraster(j).times, i*ones(size(data.tastant_2.psth_raster.spikeraster(j).times)),6,'mv','filled');
%             hold on;
%         end
%         for j=1:length(data.tastant_3.psth_raster.spikeraster)
%             i=i+1;
%             h3=scatter(data.tastant_3.psth_raster.spikeraster(j).times, i*ones(size(data.tastant_3.psth_raster.spikeraster(j).times)),6,'bv','filled');
%             hold on
%         end
%         legend([h1,h2,h3],{x{1},x{2},x{3}})
    case 4
        tastant_1                                                =timestamps1(data1==a(1));
        data.tastant_1.dig                                       =tastant_1(1:2:length(tastant_1));
%         [data.tastant_1.psth_raster] = spike2eventRasteandPSTH_NP (data.centSp.dig, data.tastant_1.dig, 100, -5000, 5000);
        tastant_2                                                =timestamps1(data1==a(2));
        data.tastant_2.dig                                       =tastant_2(1:2:length(tastant_2 ));
%         [data.tastant_2.psth_raster] = spike2eventRasteandPSTH_NP (data.centSp.dig, data.tastant_2.dig, 100, -5000, 5000);
        tastant_3                                                =timestamps1(data1==a(3));
        data.tastant_3.dig                                       =tastant_3(1:2:length(tastant_3 ));
%         [data.tastant_3.psth_raster] = spike2eventRasteandPSTH_NP (data.centSp.dig, data.tastant_3.dig, 100, -5000, 5000);
        tastant_4                                                =timestamps1(data1==a(4));
        data.tastant_4.dig                                       =tastant_4(1:2:length(tastant_4));
%         [data.tastant_4.psth_raster] = spike2eventRasteandPSTH_NP (data.centSp.dig, data.tastant_4.dig, 100, -5000, 5000);
        data.tastant_1.id=x{1};
        data.tastant_1.line=tasteID(1)+1; %get the port
        data.tastant_1.LR=y(1);
        data.leftID=[];
        data.rightID=[];
        if y(1)==1
            data.leftID='tastant_1';
        else
            data.rightID='tastant_1';
        end
        data.tastant_2.id=x{2};
        data.tastant_2.line=tasteID(2)+1; % get the port
        data.tastant_2.LR=y(2);
        if y(2)==1
            if ~isempty(data.leftID)
                data.leftID={data.leftID,'tastant_2'};
            else
                data.leftID='tastant_2';
            end
        else
            if ~isempty(data.rightID)
                data.rightID={data.rightID,'tastant_2'};
            else
                data.rightID='tastant_2';
            end
        end
        data.tastant_3.id=x{3};
        data.tastant_3.line=tasteID(3)+1; % get the port
        data.tastant_3.LR=y(3);
        if y(3)==1
            if ~isempty(data.leftID)
                data.leftID={data.leftID,'tastant_3'};
            else
                data.leftID='tastant_3';
            end
        else
            if ~isempty(data.rightID)
                data.rightID={data.rightID,'tastant_3'};
            else
                data.rightID='tastant_3';
            end
        end
        data.tastant_4.id=x{4};
        data.tastant_4.LR=y(4);
        data.tastant_4.line=tasteID(4)+1; % get the port 
        if y(4)==1
            if ~isempty(data.leftID)
                data.leftID={data.leftID,'tastant_4'};
            else
                data.leftID='tastant_4';
            end
        else
            if ~isempty(data.rightID)
                data.rightID={data.rightID,'tastant_4'};
            else
                data.rightID='tastant_4';
            end
        end
        
%         for i=1:length(data.tastant_1.psth_raster.spikeraster)
%             h1=scatter(data.tastant_1.psth_raster.spikeraster(i).times, i*ones(size(data.tastant_1.psth_raster.spikeraster(i).times)),6,'cv','filled');
%             hold on
%         end
%         for j=1:length(data.tastant_2.psth_raster.spikeraster)
%             i=i+1;
%             h2=scatter(data.tastant_2.psth_raster.spikeraster(j).times, i*ones(size(data.tastant_2.psth_raster.spikeraster(j).times)),6,'mv','filled');
%             hold on
%         end
%         for j=1:length(data.tastant_3.psth_raster.spikeraster)
%             i=i+1;
%             h3=scatter(data.tastant_3.psth_raster.spikeraster(j).times, i*ones(size(data.tastant_3.psth_raster.spikeraster(j).times)),6,'bv','filled');
%             hold on
%         end
%         for j=1:length(data.tastant_4.psth_raster.spikeraster)
%             i=i+1;
%             h4=scatter(data.tastant_4.psth_raster.spikeraster(j).times, i*ones(size(data.tastant_4.psth_raster.spikeraster(j).times)),6,'rv','filled');
%             hold on
%         end
%         legend([h1,h2,h3,h4],{x{1},x{2},x{3},x{4}})
    otherwise
        error('You have put more than 4 tastant. You need to modify the code to process it')
end

save('data.mat','data')        
performance_TAC_v3(x,tasteID)   