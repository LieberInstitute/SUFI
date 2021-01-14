%% This is a demo script to setup variables and run unmixing of lambda stack images.

clc
clear all
close   all
rng(5,'twister')

%% Load toolbox and data
toolbox = 'toolbox'; % path to toolbox
addpath(genpath(toolbox))

filename = 'example_data/BR1531_63X_MBP520_SLC17A7570_GAD1620_SNAP25690_3.czi'; % path to multiplex lambda stack that needs to be unmixed
channel = {'DAPI', 'Opal520_Lp20', 'Opal570_Lp10', 'Opal620_Lp10', 'Opal690_Lp30', 'Lipofuscin'}; % channels present in the multiplex lambda stack
singlepos = 'example_data/singlepos.csv'; % csv file with channel name and corresponding single positive .czi files - see readme for example

% Initialize pool to use multiple processor cores on the computer
npool = feature('numcores');

clus = gcp('nocreate'); % If no pool, do not create new one.
if isempty(clus)
    c = parcluster('local');
    c.NumWorkers = npool;
    parpool(c, c.NumWorkers);
end

%% FCLSU
tic
[M0, unmixed_FCLSU, stack, lambda] = main_FCLSU(filename, singlepos);
time_FCLSU = toc;
% Save results
save([filename(1:end-4),'_fclsu_unmixing.mat'],'unmixed_FCLSU', 'M0', 'lambda', 'time_FCLSU')
% Plot results
main_figure(unmixed_FCLSU, M0, lambda, channel, filename)

%{
%% ELMM - can be commented out if you prefer to use only FCLSU
tic
unmixed_ELMM = main_ELMM(unmixed_FCLSU, M0, stack);
time_ELMM = toc;
save([filename(1:end-4),'_elmm_unmixing.mat'],'unmixed_ELMM', 'M0', 'lambda', 'time_ELMM')

%% GELMM - can be commented out if you prefer to use only FCLSU
tic
unmixed_GELMM = main_GELMM(unmixed_FCLSU, M0, stack);
time_GELMM = toc;
save([filename(1:end-4),'_gelmm_unmixing.mat'],'unmixed_GELMM', 'M0', 'lambda', 'time_GELMM')
%}