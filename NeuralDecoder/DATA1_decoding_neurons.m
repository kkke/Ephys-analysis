% Decode single units, find selectivity and decoding latencies
% Plot results
%
% Luca Mazzucato August 2017
clear all;
cd('/Users/robertovincis/Documents/FSU/MatlabCodes/Data'); % change this folder.

% PREAMBLE
% this loads structure 'spike', containing 1 single neuron.
% structure: spike.(trig)(trial).spk=spike times aligned to stimulus
% delivery in each trial for stimulus (trig)
%filespikes='spike_example2.mat';
%data=load(filespikes);
%% This is how you generate the "spikes" file that will be used by the decoder:
% lest say that S, N, C, Q, and W are the "psth" matricies where
% information about the spike raster for Sucrose, NaCl, Citric Acid and
% Quinine is stored. You generate "spikes" as follow:
% spikes.W = W.spikeraster';
spikes.S = tastepsth_bin25(1).S.spikeraster';
spikes.M = tastepsth_bin25(1).M.spikeraster';
spikes.Q = tastepsth_bin25(1).Q.spikeraster';
spikes.O = tastepsth_bin25(1).O.spikeraster';
% spikes.S = tastepsth_bin25(1).S.scmatrix;
% spikes.M = tastepsth_bin25(1).M.scmatrix;
% spikes.Q = tastepsth_bin25(1).Q.scmatrix;
% spikes.O = tastepsth_bin25(1).O.scmatrix;
% 'spikeS' has fields: 'C', 'W', 'Q', 'S', 'N' -> trig
%%
trig=fieldnames(spikes);
nev=numel(trig);
%%
% decoding parameters
kind='units';
xval=2; % number of leave-out trials for each x-validation run
binsize=0.05; % size of decoding window
% windows=[0 2.5]; % total trial interval to be decoded aligned at t=0.
windows=[0 2.5]; % total trial interval to be decoded aligned at t=0.
nboot=1000; % number of bootstrap runs for each shuffled stimulus
filesave_dec='decoding'; % prefix of files to be saved
%%
[results, options] = fun_decode_units_selectivity;

%%
%--------
% PLOT
%--------
numfig=1;
fun_decode_units_selectivity_plot;

