clc
clear all
close all

toolbox = '/neural_plasticity/Vijay/Unmixing/DLPFC/Scripts/toolbox';
addpath(genpath(toolbox))

method = {'fclsu'; 'elmm'; 'gelmm'};
channel = {'DAPI', 'Opal520_Lp20', 'Opal570Lp1_0', 'Opal620_Lp10', 'Opal690_Lp30', 'Lipofuscin'};
% Load FCLSU, ELMM, GELMM and Zen unmixed files
filename = '/neural_plasticity/Vijay/Unmixing/DLPFC/BR1531_63X_MBP520_SLC17A7570_GAD1620_SNAP25690_3.czi';

% Dice similarity
segm_ZEN = load([filename(1:end-4),'_linear_unmixing_segmentation.mat'], 'Segmentations');
dice_stack = [];
for i=1:length(method)
    disp(i)
    segm = load([filename(1:end-4),'_',method{i},'_unmixing_segmentation.mat'], 'Segmentations');
    for j=1:size(channel,2)
        dice_stack(i,j) = dice(segm.Segmentations.(channel{j}), segm_ZEN.Segmentations.(channel{j}));
    end
end

% SSIM & RMSE
img_ZEN = load([filename(1:end-4),'_linear_unmixing_img.mat'], 'img');
ssim_stack = [];
rmse_stack = [];
for i=1:length(method)
    disp(i)
    img = load([filename(1:end-4),'_',method{i},'_unmixing_img.mat'], 'img');
    for j=1:size(channel,2)
        CH1 = img.img.(channel{j});
        CH2 = img_ZEN.img.(channel{j})/65535;
        ssim_stack(i,j) = ssim(CH1, CH2);
        rmse_stack(i,j) = sqrt(mean((CH1(:)- CH2(:)).^2));
    end
end

%time_stack = [time_FCLSU/60; time_ELMM/60; time_GELMM/60];
%writematrix(time_stack, [filename(1:end-4),'_time.csv'])
writematrix(rmse_stack, [filename(1:end-4),'_rmse.csv'])
writematrix(ssim_stack, [filename(1:end-4),'_ssim.csv'])
writematrix(dice_stack, [filename(1:end-4),'_dice.csv'])
