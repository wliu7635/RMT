function [imgDenoised, rank, sigma, sigmaVST] = ...
    fmriDenoise_SSVD2( ...
        epi_noisy, ...
        ks, ...
        ws, ...
        VST_ABC, ...
        wantGaussWeighting, ...
        wantRankWeighting, ...
        wantVST, ...
        k0)
% FMRIDENOISE_SSVD2
%
% Wrapper for magnitude-domain fMRI RMT denoising with optional
% Rice variance-stabilizing transformation.
%
% The final tested workflow performs:
%
%   1. Spatial Rician-noise estimation
%   2. Forward Rice variance-stabilizing transformation
%   3. Patch-based OP-SSVD denoising
%   4. Exact unbiased inverse Rice VST
%
% The last dimension of the input image is interpreted as the repeated
% measurement or time dimension.
%
%
% USAGE
% -----
%
% [imgDenoised, rank, sigma, sigmaVST] = ...
%     fmriDenoise_SSVD2( ...
%         epi_noisy, ...
%         ks, ...
%         ws, ...
%         VST_ABC, ...
%         wantGaussWeighting, ...
%         wantRankWeighting, ...
%         wantVST, ...
%         k0);
%
%
% INPUTS
% ------
%
% epi_noisy:
%   Noisy magnitude-domain EPI image series.
%
%   Expected dimensions:
%
%       [x, y, z, time]
%
%
% ks:
%   Spatial kernel used for:
%
%       estimate_noise_vst3
%       perform_riceVST3
%       perform_riceVST_EUI3
%
%   Tested value:
%
%       [7 7 3]
%
%
% ws:
%   Spatial OP-SSVD patch window.
%
%   Tested value:
%
%       [7 7 3]
%
%
% VST_ABC:
%   Rice variance-stabilizing transform used for the forward and
%   inverse VST operations.
%
%   Supported values depend on the supplied Rice-VST implementation.
%
%   Tested value:
%
%       'A'
%
%   IMPORTANT:
%
%   The successful workflow documented in this repository uses:
%
%       Noise-estimation VST method:  'B'
%       Forward VST method:           'A'
%       Inverse VST method:           'A'
%
%   The noise-estimation method is explicitly set to 'B' below to
%   preserve the tested workflow.
%
%
% wantGaussWeighting:
%   Retained for compatibility with the original wrapper interface.
%
%   Tested value:
%
%       0
%
%   IMPORTANT:
%
%   In the current OP implementation of denoise_ssvd.m, Gaussian
%   weighting is disabled. The OP overlap aggregation uses the
%   whole-patch weight:
%
%       1 / (1 + estimated rank)
%
%
% wantRankWeighting:
%   Retained for compatibility with the original wrapper interface.
%
%   Tested value:
%
%       1
%
%   IMPORTANT:
%
%   In the current OP implementation of denoise_ssvd.m, whole-patch
%   rank weighting is always applied as:
%
%       1 / (1 + estimated rank)
%
%
% wantVST:
%   Whether to apply Rician correction using forward and exact unbiased
%   inverse Rice VST.
%
%   Tested value:
%
%       1
%
%
% k0:
%   Moment orders used by the SSVD multiple-criteria rank/noise
%   estimator.
%
%   Tested value:
%
%       1:10
%
%
% OUTPUTS
% -------
%
% imgDenoised:
%   Final denoised 4D magnitude-domain image.
%
%   Dimensions:
%
%       [x, y, z, time]
%
%
% rank:
%   Patch-level rank estimates returned by denoise_ssvd.m.
%
%   Dimensions:
%
%       [x, y, z]
%
%   Values are stored at patch starting coordinates in the current
%   implementation.
%
%
% sigma:
%   Patch-level OP-SSVD noise estimates returned by denoise_ssvd.m.
%
%   Dimensions:
%
%       [x, y, z]
%
%   Values are stored at patch starting coordinates in the current
%   implementation.
%
%
% sigmaVST:
%   Spatial Rician-noise map used by the forward and inverse VST.
%
%   Dimensions:
%
%       [x, y, z]
%
%
% SUCCESSFULLY TESTED CONFIGURATION
% ---------------------------------
%
% The calling file fmriDenoising_demo.m used:
%
%   ks                  = [7 7 3]
%   ws                  = [7 7 3]
%   VST_ABC             = 'A'
%   wantGaussWeighting  = 0
%   wantRankWeighting   = 1
%   wantVST             = 1
%   k0                  = 1:10
%
% The OP-SSVD stage used:
%
%   Spatial patch:             [7 7 3]
%   Patch step:                1
%   Rank/noise estimator:      SSVD multiple criteria
%   Singular-value shrinkage:  operator-norm optimal shrinkage
%   Patch aggregation:         1 / (1 + estimated rank)
%
%
% ECHO-SPECIFIC EXECUTION
% -----------------------
%
% This tested workflow processes one echo at a time.
%
% Before running Echo 1, Echo 2, or Echo 3, manually replace:
%
%   ① cacheDir
%   ② echoLabel
%
% The echo label must correspond to the input echo selected in:
%
%   fmriDenoising_demo.m
%
% Examples:
%
% Echo 1:
%
%   echoLabel = 'echo1';
%
% Echo 2:
%
%   echoLabel = 'echo2';
%
% Echo 3:
%
%   echoLabel = 'echo3';
%
% Do not reuse one echo's cache files for another echo.
%
% If any RMT parameter is changed, delete the corresponding cache files
% before rerunning that echo. Otherwise, this function will load the
% previously completed result from the cache.
%
%
% REQUIRED FUNCTIONS
% ------------------
%
%   estimate_noise_vst3.m
%   perform_riceVST3.m
%   denoise_ssvd.m
%   perform_riceVST_EUI3.m
%   riceVST_sigmaEst.m
%   riceVST.m
%   riceVST_EUI.m
%   ssvd.m
%   squish.m
%
%
% SEE ALSO
% --------
%
%   fmriDenoising_demo
%   estimate_noise_vst3
%   perform_riceVST3
%   denoise_ssvd
%   ssvd
%   perform_riceVST_EUI3
%
%
% ORIGINAL CODE ATTRIBUTION
% -------------------------
%
% Copyright (C) 2019 CMRR at UMN
%
% Original authors:
%
%   Xiaoping Wu
%   Wei Zhu
%
% Original creation date:
%
%   September 3, 2019
%
%
% WORKFLOW MODIFICATIONS DOCUMENTED IN THIS REPOSITORY
% ---------------------------------------------------
%
%   Low-memory serial noise estimation
%   Low-memory serial forward Rice VST
%   Low-memory direct OP-SSVD patch accumulation
%   Operator-norm singular-value shrinkage integration
%   Whole-patch 1/(1+rank) overlap weighting
%   Low-memory exact unbiased inverse Rice VST
%   Echo-specific caching
%
% =====================================================================


%% Validate the input image

if nargin < 1 || isempty(epi_noisy)
    error('epi_noisy must be provided.');
end

if ndims(epi_noisy) ~= 4
    error( ...
        ['epi_noisy must be a 4D array with dimensions ', ...
         '[x, y, z, time].']);
end

if any(~isfinite(double(epi_noisy(:))))
    error('epi_noisy contains NaN or Inf values.');
end

imgSize = size(epi_noisy);


%% Set default processing parameters

if ~exist('ks', 'var') || isempty(ks)

    ks = [ ...
        min(15, imgSize(1)), ...
        min(15, imgSize(2)), ...
        min(15, imgSize(3)) ...
    ];
end

if ~exist('ws', 'var') || isempty(ws)

    ws = [3 3 3];
end

if ~exist('VST_ABC', 'var') || isempty(VST_ABC)

    VST_ABC = 'B';
end

if ~exist('wantGaussWeighting', 'var') || ...
        isempty(wantGaussWeighting)

    wantGaussWeighting = 1;
end

if ~exist('wantRankWeighting', 'var') || ...
        isempty(wantRankWeighting)

    wantRankWeighting = 0;
end

if ~exist('wantVST', 'var') || isempty(wantVST)

    wantVST = 1;
end

if ~exist('k0', 'var') || isempty(k0)

    k0 = 2;
end


%% Validate processing parameters

if isscalar(ks)
    ks = repmat(ks, 1, 3);
end

if isscalar(ws)
    ws = repmat(ws, 1, 3);
end

if numel(ks) ~= 3
    error('ks must be a scalar or a three-element vector [kx ky kz].');
end

if numel(ws) ~= 3
    error('ws must be a scalar or a three-element vector [wx wy wz].');
end

ks = double(ks(:).');
ws = double(ws(:).');

if any(ks < 1) || any(mod(ks, 1) ~= 0)
    error('All ks dimensions must be positive integers.');
end

if any(ws < 1) || any(mod(ws, 1) ~= 0)
    error('All ws dimensions must be positive integers.');
end

if any(ks > imgSize(1:3))
    error( ...
        ['The VST kernel [%d %d %d] exceeds the image dimensions ', ...
         '[%d %d %d].'], ...
        ks(1), ...
        ks(2), ...
        ks(3), ...
        imgSize(1), ...
        imgSize(2), ...
        imgSize(3));
end

if any(ws > imgSize(1:3))
    error( ...
        ['The OP-SSVD patch [%d %d %d] exceeds the image ', ...
         'dimensions [%d %d %d].'], ...
        ws(1), ...
        ws(2), ...
        ws(3), ...
        imgSize(1), ...
        imgSize(2), ...
        imgSize(3));
end


%% Display the active configuration

fprintf('\n');
fprintf('============================================================\n');
fprintf('fmriDenoise_SSVD2\n');
fprintf('============================================================\n');

fprintf( ...
    'Input dimensions: %d x %d x %d x %d\n', ...
    imgSize(1), ...
    imgSize(2), ...
    imgSize(3), ...
    imgSize(4));

fprintf( ...
    'Noise/VST kernel: [%d %d %d]\n', ...
    ks(1), ...
    ks(2), ...
    ks(3));

fprintf( ...
    'OP-SSVD patch: [%d %d %d]\n', ...
    ws(1), ...
    ws(2), ...
    ws(3));

fprintf('Forward/inverse VST type: %s\n', VST_ABC);
fprintf('VST enabled: %d\n', logical(wantVST));
fprintf('Gaussian weighting request: %s\n', ...
    mat2str(wantGaussWeighting));
fprintf('Rank weighting request: %s\n', ...
    mat2str(wantRankWeighting));
fprintf('SSVD moment orders: %s\n', mat2str(k0));

fprintf( ...
    ['OP aggregation note: denoise_ssvd.m uses whole-patch ', ...
     '1/(1+rank) weighting.\n']);

fprintf('============================================================\n\n');


%% ---------------------------------------------------------------------
% Echo-specific cache configuration
%
% Replace:
%
%   ① cacheDir
%   ② echoLabel
%
% The tested final code snapshot used Echo 3:
%
%   cacheDir = 'D:\RMT\test result';
%   echoLabel = 'echo3';
%
% To run Echo 1:
%
%   echoLabel = 'echo1';
%
% To run Echo 2:
%
%   echoLabel = 'echo2';
%
% To run Echo 3:
%
%   echoLabel = 'echo3';
%
% The echo label must match the input echo in fmriDenoising_demo.m.
%
% These two lines are the only echo-specific settings in this function.
%% ---------------------------------------------------------------------


% ① Cache/output directory
cacheDir = 'D:\RMT\test result';


% ② Echo label
%
% Replace with:
%
%   'echo1'
%   'echo2'
%   'echo3'
%
% The current GitHub example shows the final tested Echo 3 setup.

echoLabel = 'echo3';


%% Create the cache directory if needed

if ~exist(cacheDir, 'dir')

    mkdir(cacheDir);
end


%% Construct the three echo-specific cache filenames

sigmaCacheFile = fullfile( ...
    cacheDir, ...
    [echoLabel '_sigmaVST_cache.mat']);

vstCacheFile = fullfile( ...
    cacheDir, ...
    [echoLabel '_forwardVST_cache.mat']);

ssvdCacheFile = fullfile( ...
    cacheDir, ...
    [echoLabel '_OPSSVD_cache.mat']);


%% Display the active cache configuration

fprintf('Echo-specific cache configuration:\n');
fprintf('  Echo label: %s\n', echoLabel);
fprintf('  Sigma-VST cache:\n    %s\n', sigmaCacheFile);
fprintf('  Forward-VST cache:\n    %s\n', vstCacheFile);
fprintf('  OP-SSVD cache:\n    %s\n\n', ssvdCacheFile);


%% ---------------------------------------------------------------------
% Stage 1: Spatial Rician-noise estimation and forward Rice VST
%% ---------------------------------------------------------------------

if wantVST

    %% Load or estimate the spatial Rician-noise map

    if exist(sigmaCacheFile, 'file')

        fprintf( ...
            'Loading cached VST noise map:\n%s\n', ...
            sigmaCacheFile);

        cachedSigma = load( ...
            sigmaCacheFile, ...
            'sigmaVST');

        if ~isfield(cachedSigma, 'sigmaVST')
            error( ...
                ['The sigma cache does not contain the variable ', ...
                 'sigmaVST:\n%s'], ...
                sigmaCacheFile);
        end

        sigmaVST = cachedSigma.sigmaVST;

        if ~isequal(size(sigmaVST), imgSize(1:3))
            error( ...
                ['Cached sigmaVST dimensions do not match the ', ...
                 'current echo. Delete the cache and rerun:\n%s'], ...
                sigmaCacheFile);
        end

    else

        fprintf( ...
            ['No cached VST noise map was found.\n', ...
             'Estimating the spatial Rician-noise map...\n']);

        % IMPORTANT:
        %
        % The successful tested workflow used VST estimator B for
        % spatial noise estimation, while the forward and inverse
        % transforms used the caller-supplied VST_ABC value.
        %
        % With the tested fmriDenoising_demo.m configuration:
        %
        %   Noise estimation:  B
        %   Forward VST:       A
        %   Inverse VST:       A

        sigmaVST = estimate_noise_vst3( ...
            epi_noisy, ...
            ks, ...
            'B');

        if any(~isfinite(double(sigmaVST(:))))
            error( ...
                'The estimated sigmaVST map contains NaN or Inf.');
        end

        save( ...
            sigmaCacheFile, ...
            'sigmaVST', ...
            '-v7.3');

        fprintf( ...
            'Saved VST noise-map cache:\n%s\n', ...
            sigmaCacheFile);
    end


    %% Load or calculate the forward Rice VST output

    if exist(vstCacheFile, 'file')

        fprintf( ...
            'Loading cached forward-VST data:\n%s\n', ...
            vstCacheFile);

        cachedVST = load( ...
            vstCacheFile, ...
            'imgRaw');

        if ~isfield(cachedVST, 'imgRaw')
            error( ...
                ['The forward-VST cache does not contain the ', ...
                 'variable imgRaw:\n%s'], ...
                vstCacheFile);
        end

        imgRaw = cachedVST.imgRaw;

        if ~isequal(size(imgRaw), imgSize)
            error( ...
                ['Cached forward-VST dimensions do not match the ', ...
                 'current echo. Delete the cache and rerun:\n%s'], ...
                vstCacheFile);
        end

    else

        fprintf( ...
            ['No forward-VST cache was found.\n', ...
             'Running the forward Rice VST...\n']);

        imgRaw = perform_riceVST3( ...
            epi_noisy, ...
            sigmaVST, ...
            ks, ...
            VST_ABC);

        if ~isequal(size(imgRaw), imgSize)
            error( ...
                ['perform_riceVST3 returned an unexpected ', ...
                 'output size.']);
        end

        if any(~isfinite(double(imgRaw(:))))
            error( ...
                'The forward-VST output contains NaN or Inf.');
        end

        save( ...
            vstCacheFile, ...
            'imgRaw', ...
            '-v7.3');

        fprintf( ...
            'Saved forward-VST cache:\n%s\n', ...
            vstCacheFile);
    end

else

    %% VST-disabled branch

    fprintf( ...
        ['VST is disabled. OP-SSVD will operate directly on the ', ...
         'input image series.\n']);

    imgRaw = epi_noisy;
    sigmaVST = [];
end


%% ---------------------------------------------------------------------
% Stage 2: OP-SSVD RMT denoising
%
% The tested workflow uses:
%
%   Patch window:             ws = [7 7 3]
%   Spatial step:            1
%   Rank/noise estimator:    SSVD multiple criteria
%   Moment orders:           k0 = 1:10
%   Shrinkage:               operator-norm optimal shrinkage
%   Patch aggregation:       1 / (1 + estimated rank)
%% ---------------------------------------------------------------------

step = 1;


%% Load or calculate OP-SSVD results

if exist(ssvdCacheFile, 'file')

    fprintf( ...
        'Loading cached OP-SSVD results:\n%s\n', ...
        ssvdCacheFile);

    cachedSSVD = load( ...
        ssvdCacheFile, ...
        'imgDenoised', ...
        'rank', ...
        'sigma');

    requiredFields = { ...
        'imgDenoised', ...
        'rank', ...
        'sigma' ...
    };

    for fieldIndex = 1:numel(requiredFields)

        if ~isfield( ...
                cachedSSVD, ...
                requiredFields{fieldIndex})

            error( ...
                ['The OP-SSVD cache is missing the variable ', ...
                 '%s:\n%s'], ...
                requiredFields{fieldIndex}, ...
                ssvdCacheFile);
        end
    end

    imgDenoised = cachedSSVD.imgDenoised;
    rank = cachedSSVD.rank;
    sigma = cachedSSVD.sigma;

    if ~isequal(size(imgDenoised), imgSize)
        error( ...
            ['Cached OP-SSVD output dimensions do not match the ', ...
             'current echo. Delete the cache and rerun:\n%s'], ...
            ssvdCacheFile);
    end

else

    fprintf( ...
        ['No OP-SSVD cache was found.\n', ...
         'Running low-memory OP-SSVD denoising...\n']);

    [imgDenoised, rank, sigma] = denoise_ssvd( ...
        imgRaw, ...
        ws, ...
        step, ...
        wantGaussWeighting, ...
        wantRankWeighting, ...
        k0);

    if ~isequal(size(imgDenoised), imgSize)
        error( ...
            ['denoise_ssvd returned an unexpected output ', ...
             'size.']);
    end

    if any(~isfinite(double(imgDenoised(:))))
        error( ...
            'The OP-SSVD output contains NaN or Inf.');
    end

    save( ...
        ssvdCacheFile, ...
        'imgDenoised', ...
        'rank', ...
        'sigma', ...
        '-v7.3');

    fprintf( ...
        'Saved OP-SSVD cache:\n%s\n', ...
        ssvdCacheFile);
end


%% ---------------------------------------------------------------------
% Stage 3: Exact unbiased inverse Rice VST
%
% This stage converts the OP-SSVD result from the variance-stabilized
% domain back to magnitude-domain image intensity.
%% ---------------------------------------------------------------------

if wantVST

    fprintf( ...
        'Running the exact unbiased inverse Rice VST...\n');

    imgDenoised = perform_riceVST_EUI3( ...
        imgDenoised, ...
        sigmaVST, ...
        ks, ...
        VST_ABC);

    if ~isequal(size(imgDenoised), imgSize)
        error( ...
            ['perform_riceVST_EUI3 returned an unexpected ', ...
             'output size.']);
    end

    if any(~isfinite(double(imgDenoised(:))))
        error( ...
            ['The final inverse-VST output contains ', ...
             'NaN or Inf.']);
    end
end


%% ---------------------------------------------------------------------
% Final validation and reporting
%% ---------------------------------------------------------------------

if any(~isfinite(double(imgDenoised(:))))
    error( ...
        'The final denoised output contains NaN or Inf.');
end

if wantVST && min(imgDenoised(:)) < 0
    error( ...
        ['The final magnitude-domain output contains negative ', ...
         'values after inverse VST.']);
end

fprintf('\n');
fprintf('============================================================\n');
fprintf('fmriDenoise_SSVD2 completed successfully\n');
fprintf('============================================================\n');

fprintf('Echo label: %s\n', echoLabel);

fprintf( ...
    'Final dimensions: %d x %d x %d x %d\n', ...
    size(imgDenoised, 1), ...
    size(imgDenoised, 2), ...
    size(imgDenoised, 3), ...
    size(imgDenoised, 4));

fprintf( ...
    'Final intensity range: minimum=%g, maximum=%g\n', ...
    double(min(imgDenoised(:))), ...
    double(max(imgDenoised(:))));

if ~isempty(rank)

    fprintf( ...
        'Patch-rank range: minimum=%g, maximum=%g\n', ...
        double(min(rank(:))), ...
        double(max(rank(:))));
end

if ~isempty(sigma)

    fprintf( ...
        'OP-SSVD sigma range: minimum=%g, maximum=%g\n', ...
        double(min(sigma(:))), ...
        double(max(sigma(:))));
end

if ~isempty(sigmaVST)

    fprintf( ...
        'VST sigma range: minimum=%g, maximum=%g\n', ...
        double(min(sigmaVST(:))), ...
        double(max(sigmaVST(:))));
end

fprintf('============================================================\n\n');

end
