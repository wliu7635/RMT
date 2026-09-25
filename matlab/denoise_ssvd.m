function [im_out, rank_out, noise_out] = denoise_ssvd( ...
    nim, b, step, wantGaussWeighting, wantRankWeighting, k0)
% DENOISE_SSVD
% Low-memory patch-based fMRI RMT denoising using:
%
%   1. The 2022 fMRI SSVD multiple-criteria rank/noise estimator
%   2. Operator-norm optimal singular-value shrinkage ('op')
%   3. The 2020 whole-patch rank weighting:
%
%          weight = 1 / (1 + estimated rank)
%
% This implementation processes one patch at a time and directly
% accumulates the reconstructed patch estimates. It avoids the original
% large five-dimensional PARFOR arrays.
%
% IMPORTANT
% ---------
% wantGaussWeighting and wantRankWeighting are retained in the function
% signature for compatibility with fmriDenoise_SSVD2.m.
%
% For this OP implementation, the aggregation follows the supplied 2020
% denoise_optim_SVShrinkage3.m implementation:
%
%       wei = 1/(1+r);
%       W   = wei*ones(size(X));
%       X   = X*wei;
%
% Therefore Gaussian delta weighting is not used in this function.
%
% Inputs
% ------
% nim:
%   Input VST-domain 4D data [x y z M].
%
% b:
%   Spatial SSVD patch dimensions. Scalar or [bx by bz].
%
% step:
%   Spatial step between neighboring patches.
%
% wantGaussWeighting:
%   Retained for interface compatibility. Not used in this OP aggregation.
%
% wantRankWeighting:
%   Retained for interface compatibility. The 2020 OP weighting
%   1/(1+R) is always used.
%
% k0:
%   Moment orders used by the 2022 SSVD rank/noise estimator.
%
% Outputs
% -------
% im_out:
%   Denoised VST-domain 4D data.
%
% rank_out:
%   Patch rank values stored at patch starting positions.
%
% noise_out:
%   Patch noise estimates stored at patch starting positions.
%
% Original algorithm copyright:
% Copyright (C) 2019 CMRR at UMN
% Authors: Xiaoping Wu and Wei Zhu
%
% Low-memory and OP integration modification:
% Direct serial patch accumulation, 2022 rank/noise estimation,
% operator-norm shrinkage, and 2020 rank-weighted patch aggregation.

%% Defaults

if nargin < 2 || isempty(b)
    b = 5;
end

if isscalar(b)
    b = repmat(b, 1, 3);
end

if numel(b) ~= 3
    error('b must be a scalar or a three-element vector [bx by bz].');
end

b = double(b(:).');

bx = b(1);
by = b(2);
bz = b(3);

if any(b < 1) || any(mod(b,1) ~= 0)
    error('All SSVD window dimensions must be positive integers.');
end

if nargin < 3 || isempty(step)
    step = 1;
end

if ~isscalar(step) || step < 1 || mod(step,1) ~= 0
    error('step must be a positive integer.');
end

% Retained only for interface compatibility.
if nargin < 4 || isempty(wantGaussWeighting)
    wantGaussWeighting = [0 0];
end

% Retained only for interface compatibility.
if nargin < 5 || isempty(wantRankWeighting)
    wantRankWeighting = 1;
end

if nargin < 6 || isempty(k0)
    k0 = 2;
end

%% Input preparation

[sx, sy, sz, M] = size(nim);

if bx > sx || by > sy || bz > sz
    error(['SSVD window [%d %d %d] exceeds image dimensions ', ...
           '[%d %d %d].'], ...
           bx, by, bz, sx, sy, sz);
end

time0 = clock;

fprintf('--------start OP-SSVD denoising--------\n');

fprintf('Input dimensions: %d x %d x %d x %d\n', ...
    sx, sy, sz, M);

fprintf('SSVD window: %d x %d x %d\n', ...
    bx, by, bz);

fprintf('Step: %d\n', step);

fprintf(['Shrinkage: operator-norm optimal shrinkage\n', ...
         'Rank/noise estimator: SSVD multiple criteria\n', ...
         'Patch aggregation: whole-patch 1/(1+R) weighting\n']);

if any(wantGaussWeighting)
    fprintf(['Note: wantGaussWeighting was supplied as [%s], ', ...
             'but Gaussian weighting is disabled for this OP ', ...
             'whole-patch aggregation.\n'], ...
             num2str(wantGaussWeighting));
end

if ~wantRankWeighting
    fprintf(['Note: wantRankWeighting was supplied as 0, but the ', ...
             '2020 OP aggregation always applies 1/(1+R).\n']);
end

%% Allocate direct low-memory accumulators

% One complete 4D accumulator.
% Double precision preserves weighted accumulation accuracy.
Ys = zeros(sx, sy, sz, M, 'double');

% The patch weight is identical across time, so one 3D denominator is
% mathematically equivalent to a replicated 4D denominator.
W = zeros(sx, sy, sz, 'double');

% Patch rank and noise maps.
R = zeros(sx, sy, sz, 'double');
S = zeros(sx, sy, sz, 'double');

%% Patch starting positions

indices_x = [1:step:sx-bx, sx-bx+1];
indices_y = [1:step:sy-by, sy-by+1];
indices_z = [1:step:sz-bz, sz-bz+1];

indices_x = unique(indices_x, 'stable');
indices_y = unique(indices_y, 'stable');
indices_z = unique(indices_z, 'stable');

nPatchX = length(indices_x);
nPatchY = length(indices_y);
nPatchZ = length(indices_z);

fprintf('Patch starting positions: %d x %d x %d\n', ...
    nPatchX, nPatchY, nPatchZ);

fprintf('Total spatial patches: %d\n', ...
    nPatchX * nPatchY * nPatchZ);

disp('-> OP-SSVD using direct low-memory patch accumulation...');

%% Process every spatial patch

for ii = 1:nPatchX

    i = indices_x(ii);

    fprintf('--- OP-SSVD x patch %d of %d, starting at x=%d ---\n', ...
        ii, nPatchX, i);

    for jj = 1:nPatchY

        j = indices_y(jj);

        for kk = 1:nPatchZ

            k = indices_z(kk);

            %% Extract one 4D spatial patch

            B1 = nim( ...
                i:i+bx-1, ...
                j:j+by-1, ...
                k:k+bz-1, ...
                :);

            %% Apply hybrid fMRI RMT-OP configuration

            % Rank/noise estimator:
            %   2022 SSVD multiple-criteria method
            %
            % Singular-value manipulation:
            %   Operator-norm optimal shrinkage ('op')
            %
            % Moment orders:
            %   k0

            [Ysp, Rp, Sigmap] = ssvd( ...
                double(B1), ...
                'op', ...
                'ssvd', ...
                k0);

            Ysp = double(Ysp);

            %% Validate patch outputs

            if ~isequal(size(Ysp), size(B1))
                error(['OP-SSVD output-size mismatch at patch ', ...
                       '[x=%d, y=%d, z=%d].'], ...
                       i, j, k);
            end

            if any(~isfinite(Ysp(:)))
                error(['OP-SSVD produced NaN or Inf at patch ', ...
                       '[x=%d, y=%d, z=%d].'], ...
                       i, j, k);
            end

            if ~isscalar(Rp) || ~isfinite(Rp) || Rp < 0
                error(['Invalid rank estimate at patch ', ...
                       '[x=%d, y=%d, z=%d].'], ...
                       i, j, k);
            end

            if ~isscalar(Sigmap) || ...
                    ~isfinite(Sigmap) || Sigmap < 0

                error(['Invalid noise estimate at patch ', ...
                       '[x=%d, y=%d, z=%d].'], ...
                       i, j, k);
            end

            %% Apply the 2020 OP whole-patch weighting

            patchScalarWeight = 1 / (1 + Rp);

            if ~isfinite(patchScalarWeight) || ...
                    patchScalarWeight <= 0

                error(['Invalid patch weight at patch ', ...
                       '[x=%d, y=%d, z=%d].'], ...
                       i, j, k);
            end

            % Every voxel in the patch receives the same positive weight.
            patchWeight = ...
                patchScalarWeight .* ...
                ones(bx, by, bz, 'double');

            weightedPatch = Ysp .* patchScalarWeight;

            %% Accumulate the weighted reconstructed patch

            Ys( ...
                i:i+bx-1, ...
                j:j+by-1, ...
                k:k+bz-1, ...
                :) = ...
            Ys( ...
                i:i+bx-1, ...
                j:j+by-1, ...
                k:k+bz-1, ...
                :) + weightedPatch;

            %% Accumulate the spatial denominator

            W( ...
                i:i+bx-1, ...
                j:j+by-1, ...
                k:k+bz-1) = ...
            W( ...
                i:i+bx-1, ...
                j:j+by-1, ...
                k:k+bz-1) + patchWeight;

            %% Store patch-level rank and noise estimates

            R(i,j,k) = Rp;
            S(i,j,k) = Sigmap;

        end
    end
end

%% Validate weighted aggregation

disp('-> validate OP-SSVD accumulated weights...');

zeroWeightMask = (W == 0);
nZeroWeight = nnz(zeroWeightMask);

fprintf('Voxels with zero accumulated weight: %d\n', ...
    nZeroWeight);

positiveWeights = W(W > 0);

if isempty(positiveWeights)
    error('All accumulated OP-SSVD weights are zero.');
end

fprintf('Positive accumulated weight minimum: %.17g\n', ...
    min(positiveWeights));

fprintf('Positive accumulated weight maximum: %.17g\n', ...
    max(positiveWeights));

% With whole-patch positive weighting and the original border-complete
% patch coordinates, every spatial voxel should have W > 0.
if nZeroWeight > 0

    diagnosticFile = ...
        'D:\RMT\test result\echo1_OPSSVD_zeroWeight_diagnostic.mat';

    save( ...
        diagnosticFile, ...
        'zeroWeightMask', ...
        'W', ...
        'indices_x', ...
        'indices_y', ...
        'indices_z', ...
        'b', ...
        'step', ...
        '-v7.3');

    error(['Unexpected zero accumulated weight for %d voxels. ', ...
           'No fallback was applied. Diagnostic data were saved to: %s'], ...
           nZeroWeight, diagnosticFile);

end

%% Exact weighted overlap aggregation

% This is mathematically equivalent to the original:
%
%     im_out = Ys ./ W
%
% where the same 3D W is replicated across all time points.
%
% Normalize one time point at a time to reduce peak memory.

for tt = 1:M

    Ys(:,:,:,tt) = ...
        Ys(:,:,:,tt) ./ W;

end

im_out = Ys;
rank_out = R;
noise_out = S;

%% Final validation and reporting

if any(~isfinite(im_out(:)))
    error('The final OP-SSVD output contains NaN or Inf.');
end

fprintf('OP-SSVD output range: minimum=%g, maximum=%g\n', ...
    min(im_out(:)), max(im_out(:)));

fprintf('Estimated patch-rank range: minimum=%g, maximum=%g\n', ...
    min(R(:)), max(R(:)));

fprintf('Estimated SSVD sigma range: minimum=%g, maximum=%g\n', ...
    min(S(:)), max(S(:)));

fprintf('Total elapsed time = %f min\n\n', ...
    etime(clock,time0) / 60);

end
