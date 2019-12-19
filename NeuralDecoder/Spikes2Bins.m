% function temp_rate=Spikes2Bins(Spikes,bins);
% INPUT Spikes(ntrials,nunits).spk, row vector of spike times
% OUTPUT temp_rate, array of rates, 1st=trials, 
%               if numel(bins)=1, 2nd=units, 
%               if numel(bins)>1 2nd=bins, 3rd=units
%
% Luca Mazzucato 
% July 2015
 
function temp_rate=Spikes2Bins(Spikes,bins)
 
gnunits=size(Spikes,2);
NBins=numel(bins)-1;
ntrials=size(Spikes,1);
BinSize=diff(bins(1:2));
if NBins>1
    temp_rate=NaN(ntrials,NBins,gnunits); % rows=bins, cols=units
elseif NBins==1
    temp_rate=NaN(ntrials,gnunits); % rows=bins, cols=units
end
for trial=1:ntrials
    spcnt_all=zeros(NBins,gnunits);
    for unit=1:gnunits
%         if ~isempty(Spikes(trial,unit).spk)
            temp_data=Spikes(trial,unit).times;
%             temp_data=Spikes(trial,unit);

            % spike counts
            spcnt=NaN(NBins,1);
            for b=1:NBins
                spcnt(b)=sum(temp_data<bins(b+1) & temp_data>=bins(b)); % spikes per bin 
            end
            spcnt_all(1:NBins,unit)=spcnt(1:NBins)/BinSize;
%         end
    end
    if NBins>1
        temp_rate(trial,1:NBins,1:gnunits)=spcnt_all(1:NBins,1:gnunits); % rows=bins, cols=units
    elseif NBins==1
        temp_rate(trial,1:gnunits)=spcnt_all(1:gnunits); % rows=bins, cols=units
    end
end