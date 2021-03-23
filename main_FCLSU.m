function [M0, unmixed, image6d, lambda] = main_FCLSU(filename, varargin)
% Function to unmix lambda stack based on Fully Constrained Linear Spectral Unmixing (FCLSU) algorithm
% Perform a Linear least squares with nonnegativity constraints.
% -------------------------------------------------------------------------
% Inputs: filename - path to lambda stack
%         varargin (singlepos) - path to a csv file with channel, single positive image location annotation
%                      example - DAPI, path to DAPI single positive
%                                Opal520, path to Opal520 single positive
% Outputs: M0 - extracted fingerprints
%          unmixed - unmixed 6D stack
%          image6d - Input image formatted to 6D stack
%          lambda - lambda sampling values of imaging (eg: 400nm, 440nm...)
% Copyright (2020) Lieber Institute for Brain Development, Baltimore MD


% check number of inputs
minnargin = 1;
maxnargin = 4;
narginchk(minnargin,maxnargin);

% check mandatory inputs

if isfile(filename)
    warning('off','all');
    out = ReadImage6D(filename);
    X = out{2}.SizeX;
    Y = out{2}.SizeY;
    Z = out{2}.SizeZ;
    S = out{2}.SeriesCount;
    T = out{2}.SizeT;
    nlambda = out{2}.SizeC;
    lambda = str2num(cell2mat(out{2}.Channels'));
    image6d = out{1}; %dims = series,time, z, c, x, y
    image6d = image6d/max(image6d(:));
    warning('on','all');
else
    error('Check input file location. No file in the directory')
end


% initialize optional inputs

singlepos = [];

% change optional parameters if provided

for i = 1:length(varargin)
    switch i
        case 1
            singlepos = varargin{i};
            singlepos = readtable(singlepos, 'ReadVariableNames', true, 'Delimiter', ',');
    end
end

%% Find fingerprints using VCA
C = size(singlepos,1) + 1;
if(size(singlepos,1))
    M0 = [];
    noise = [];
    P_vca = 2; % channel and noise
    for i = 1:size(singlepos,1)
        warning('off','all');
        img = ReadImage6D(singlepos.file{i});
        warning('on','all');
        img = img{1};
        
        r_cube = squeeze(rescale(img));
        r_cube = permute(r_cube, [2, 3, 1]);
        [m,n,L] = size(r_cube);
        r = reshape(r_cube, [m*n, L])';
        
        M_print = vca(r,'Endmembers',P_vca);
        [~,maxid]= max(mean(M_print));
        [~,minid]= min(mean(M_print));
        M0(:,i) = M_print(:, maxid);
        noise(:,i) = M_print(:, minid);
    end
    M0(:, end+1) = mean(noise, 2);
else
    P_vca = nlambda + 1; % num channels plus noise
    img = squeeze(image6d(:,:,round(Z/2),:,:,:));
    r_cube = permute(img, [2, 3, 1]);
    [m,n,L] = size(r_cube);
    r = reshape(r_cube, [m*n, L])';
    
    M0 = vca(r,'Endmembers',P_vca);
end
%% Run FCLSU
%unmixed = zeros(S, T, Z, C, X, Y);
unmixed = zeros(X, Y, C, Z);
parfor i=1:Z
    r_cube = permute(squeeze(image6d(:,:,i,:,:,:)), [2,3,1]);
    
    [m,n,L] = size(r_cube);
    r = reshape(r_cube, [m*n, L])';
    %% Fully Constrained Least Squares Unmixing (FCLSU)
    unmixed(:,:,:,i) = reshape(FCLSU(r,M0),m,n,C);
    fprintf('Done with frame: %d of %d \n', i, Z)
end
end
