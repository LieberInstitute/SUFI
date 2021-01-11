%% This is a demo script to setup variables and run unmixing of lambda stack images.

clc
clear all
close   all
rng(5,'twister')

%% Load toolbox and data
algo = "FCLSU";
toolbox = '/neural_plasticity/Vijay/Unmixing/DLPFC/Scripts/toolbox';
addpath(genpath(toolbox))

channel = {'DAPI', 'Opal520_Lp20', 'Opal570_Lp10', 'Opal620_Lp10', 'Opal690_Lp30', 'Lipofuscin'};
filename = '/neural_plasticity/Vijay/Unmixing/DLPFC/BR1531_63X_MBP520_SLC17A7570_GAD1620_SNAP25690_3.czi';
singlepos = '/neural_plasticity/Vijay/Unmixing/DLPFC/singlepos.csv';

% Initialize pool
%npool = feature('numcores');
npool = 40;

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

%% ELMM and GELMM
tic
unmixed_ELMM = main_ELMM(unmixed_FCLSU, M0, stack);
time_ELMM = toc;
save([filename(1:end-4),'_elmm_unmixing.mat'],'unmixed_ELMM', 'M0', 'lambda', 'time_ELMM')

tic
unmixed_GELMM = main_GELMM(unmixed_FCLSU, M0, stack);
time_GELMM = toc;
save([filename(1:end-4),'_gelmm_unmixing.mat'],'unmixed_GELMM', 'M0', 'lambda', 'time_GELMM')