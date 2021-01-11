function [unmixed] = main_GELMM(unmixed_FCLSU, M0, img)
% Function to unmix lambda stack based on Generalized Full Extended Linear Mixing Model (GELMM) algorithm
% -------------------------------------------------------------------------
% Inputs: unmixed_FCLSU - 6D stack of unmixed output from FCLSU
%         M0 - extracted fingerprints
%         img - Input image as 6D stack
%
% Outputs: unmixed - unmixed 6D stack using ELMM
% Copyright (2020) Lieber Institute for Brain Development, Baltimore MD

nnorm = '1,1';
verbose = true;
maxiter_anls = 20;
maxiter_admm = 100;
epsilon_s = 10^(-3);
epsilon_a = 10^(-3);
epsilon_psi = 10^(-3);
epsilon_admm_abs = 10^(-2);
epsilon_admm_rel = 10^(-2);
lambda_s   = 1;
lambda_a   = 0.01;
lambda_psi = 1e-6;
P = size(M0,2);
Z = size(unmixed_FCLSU, 4);
unmixed = zeros(size(unmixed_FCLSU));

%% Run GELMM
for i=1:Z
    r_cube = permute(squeeze(img(:,:,i,:,:,:)), [2, 3, 1]);
    [m,n,L] = size(r_cube);
    psis_init = ones(L, P, m*n);
    A_FCLSU = reshape(permute(squeeze(unmixed_FCLSU(:,:,:,i)), [3,1,2]),[P, m*n]);
    
    [A_GELMM, ~, ~, ~] = GLMM_ADMM(r_cube, A_FCLSU, psis_init, M0,lambda_s,lambda_a,lambda_psi,nnorm,verbose,maxiter_anls,maxiter_admm,epsilon_s,epsilon_a,epsilon_psi,epsilon_admm_abs,epsilon_admm_rel);
    unmixed(:,:,:,i) = reshape(A_GELMM',m,n,P);
end
end