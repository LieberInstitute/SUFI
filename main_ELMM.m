function [unmixed] = main_ELMM(unmixed_FCLSU, M0, img)
% Function to unmix lambda stack based on Full Extended Linear Mixing Model (ELMM) algorithm
% -------------------------------------------------------------------------
% Inputs: unmixed_FCLSU - 6D stack of unmixed output from FCLSU
%         M0 - extracted fingerprints
%         img - Input image as 6D stack
%
% Outputs: unmixed - unmixed 6D stack using ELMM
% Copyright (2020) Lieber Institute for Brain Development, Baltimore MD


% optional parameters
nnorm = '1,1'; % Use a Total Variation on the abundances
verbose = true; % display
maxiter_anls = 20;
maxiter_admm = 100;
epsilon_s = 10^(-3);
epsilon_a = 10^(-3);
epsilon_psi = 10^(-3);
epsilon_admm_abs = 10^(-2);
epsilon_admm_rel = 10^(-2);
lambda_s = 0.4;
lambda_a = 0.005;
lambda_psi = 0.001;
P = size(M0,2);
Z = size(unmixed_FCLSU, 4);
unmixed = zeros(size(unmixed_FCLSU));

%% Run ELMM
for i=1:Z
    r_cube = permute(squeeze(img(:,:,i,:,:,:)), [2, 3, 1]);
    [m,n,~] = size(r_cube);
    psis_init = ones(P, m*n);
    A_FCLSU = reshape(permute(squeeze(unmixed_FCLSU(:,:,:,i)), [3,1,2]),[P, m*n]);
    
    [A_ELMM, ~, ~, ~] = ELMM_ADMM(r_cube, A_FCLSU, psis_init, M0,lambda_s,lambda_a,lambda_psi,nnorm,verbose,maxiter_anls,maxiter_admm,epsilon_s,epsilon_a,epsilon_psi,epsilon_admm_abs,epsilon_admm_rel);
    unmixed(:,:,:, i) = reshape(A_ELMM',m,n,P);
end
end