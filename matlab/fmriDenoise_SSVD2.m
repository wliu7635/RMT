function [imgDenoised,rank,sigma,sigmaVST] = fmriDenoise_SSVD2(epi_noisy, ks, ws, VST_ABC, wantGaussWeighting,wantRankWeighting,wantVST,k0)

%
% 'fmriDenoise_SSVD2' is a wraper of RMT-based denoising algorithm with (or
% without) VST to denoise magnitude (or complex) 3D/4D images with the last
% dimention representing the repetitive measurement. 
%
% Usage: [imgDenoised,rank,sigma,sigmaVST] = fmriDenoise_SSVD2 (epi_noisy, ks, ws, VST_ABC, wantGaussWeighting,wantRankWeighting,wantVST,k0)
%
% Returns
% -------
% imgDenoised: [x y z M]
% rank: [x y z], ranks for individual patches after singular value
%   shrinkage.
% sigma: [x y z], estimated noise
% sigmaVST: [x,y,z], estimated noise from VST
%
% Expects
% -------
% epi_noisy: noisy epi images
% ks: kernel size used in VST operation
% ws: window size, defaults is 3x3x3 
% VST_ABC: which variance stabilizer in VST should be used, 'A' or 'B' (default)
% wantGaussWeighting: how much Gaussian weighting is used on the patch.
%   Default: [1 400]. The first number is whether to
%   use Gaussian weighting and the second number is width factor, inversely
%   proportional to the width of the window. 400 means delta weighting.
% wantRankWeighting: Whether to use rank weighting, default is 0.
%   k0: the moments for the generalized circle law to dertermine the rank
%   used in ssvd.m for 'ssvd' method.
% wantVST: whether use VST to do Rician correction
%
% See also: denoise_ssvd.m ssvd.m
%
%
% Copyright (C) 2019 CMRR at UMN
% Author: Xiaoping Wu <xpwu@cmrr.umn.edu> and Wei Zhu <zhuwei@umn.edu>
% Created: Tue Sep  3 15:02:51 2019

imgSize = size(epi_noisy);

% parameter setting
if ~exist('ks','var') || isempty(ks)
    ks = [min(15,imgSize(1)), min(15,imgSize(2)), min(15,imgSize(3))]; % kernel size for noise estimation
end
if ~exist('ws','var') || isempty(ws)
    ws = [3 3 3]; % window size for denosing
end
if ~exist('VST_ABC','var') || isempty(VST_ABC)
    VST_ABC = 'B'; % choose the function for noise estimation; VST-B recommended
end
if ~exist('wantGaussWeighting','var') || isempty(wantGaussWeighting)
  wantGaussWeighting = 1;
end
if ~exist('wantRankWeighting','var') || isempty(wantRankWeighting)
  wantRankWeighting = 0;
end
if ~exist('wantVST','var') || isempty(wantVST)
  wantVST = 1;
end
if ~exist('k0','var') || isempty(k0)
  k0 = 2;
end


% Parallel pool is disabled during low-memory noise estimation.
% Other processing functions may be converted to low-memory execution
% separately if required.

% pool = gcp('nocreate');
% if isempty(pool)
%     mypool = parpool('Processes',2);
% end

% figure; imshow(makeMontage(epi0(:,:,sliceSelect,1),1,2,'xy'),[]);
% figure; imshow(makeMontage(epi0_noisy(:,:,sliceSelect,1),1,2,'xy'),[]);


%% If VST is wanted
if wantVST

    % Cache the VST noise map so that a completed noise-estimation step
    % does not need to be repeated after a later-stage error.
    sigmaCacheFile = ...
        'D:\RMT\test result\echo3_sigmaVST_cache.mat';

    if exist(sigmaCacheFile, 'file')

        fprintf('Loading cached VST noise map:\n%s\n', ...
            sigmaCacheFile);

        cachedSigma = load(sigmaCacheFile, 'sigmaVST');
        sigmaVST = cachedSigma.sigmaVST;

    else

        fprintf('No cached VST noise map found. Estimating noise...\n');

        sigmaVST = estimate_noise_vst3( ...
            epi_noisy, ...
            ks, ...
            'B');

        if ~exist('D:\RMT\test result', 'dir')
            mkdir('D:\RMT\test result');
        end

        save( ...
            sigmaCacheFile, ...
            'sigmaVST', ...
            '-v7.3');

        fprintf('Saved VST noise map cache:\n%s\n', ...
            sigmaCacheFile);

    end

    % Forward Rice variance-stabilizing transformation.
    % Cache the forward VST output so that it does not need to be recomputed
% if SSVD or inverse VST fails later.
vstCacheFile = ...
    'D:\RMT\test result\echo3_forwardVST_cache.mat';

if exist(vstCacheFile, 'file')

    fprintf('Loading cached forward VST data:\n%s\n', ...
        vstCacheFile);

    cachedVST = load(vstCacheFile, 'imgRaw');
    imgRaw = cachedVST.imgRaw;

else

    fprintf('No forward VST cache found. Running forward VST...\n');

    imgRaw = perform_riceVST3( ...
        epi_noisy, ...
        sigmaVST, ...
        ks, ...
        VST_ABC);

    save( ...
        vstCacheFile, ...
        'imgRaw', ...
        '-v7.3');

    fprintf('Saved forward VST cache:\n%s\n', ...
        vstCacheFile);

end

else

    imgRaw = epi_noisy;
    sigmaVST = [];

end

%% denoise using SSVD

step = 1;

ssvdCacheFile = ...
    'D:\RMT\test result\echo3_OPSSVD_cache.mat';

if exist(ssvdCacheFile, 'file')

    fprintf('Loading cached SSVD results:\n%s\n', ...
        ssvdCacheFile);

    cachedSSVD = load( ...
        ssvdCacheFile, ...
        'imgDenoised', ...
        'rank', ...
        'sigma');

    imgDenoised = cachedSSVD.imgDenoised;
    rank = cachedSSVD.rank;
    sigma = cachedSSVD.sigma;

else

    fprintf('No SSVD cache found. Running SSVD denoising...\n');

    [imgDenoised, rank, sigma] = denoise_ssvd( ...
        imgRaw, ...
        ws, ...
        step, ...
        wantGaussWeighting, ...
        wantRankWeighting, ...
        k0);

    save( ...
        ssvdCacheFile, ...
        'imgDenoised', ...
        'rank', ...
        'sigma', ...
        '-v7.3');

    fprintf('Saved SSVD cache:\n%s\n', ...
        ssvdCacheFile);

end
%% If VST is used, do EUI VST
if wantVST
    % EUI VST
    imgDenoised = perform_riceVST_EUI3(imgDenoised, sigmaVST, ks, VST_ABC);    
 
end




