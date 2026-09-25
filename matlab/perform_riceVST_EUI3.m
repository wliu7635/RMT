function fz_data = perform_riceVST_EUI3(data, sigma, b, VST_ABC)
% PERFORM_RICEVST_EUI3
% Perform exact unbiased inverse variance-stabilizing transformation
% for spatially varying Rician noise.
%
% This is a low-memory serial implementation of the original
% perform_riceVST_EUI3.m.
%
% The original algorithm is preserved:
%   1. Use overlapping spatial patches.
%   2. Calculate the mean noise level inside each patch.
%   3. Apply riceVST_EUI to each patch.
%   4. Aggregate multiple patch estimates for every voxel.
%   5. Average overlapping estimates.
%
% This version avoids allocating the very large five-dimensional arrays
% data0, sigma0, Ys0, and W0.
%
% Usage:
%   fz_data = perform_riceVST_EUI3(data, sigma, b, VST_ABC)
%
% Returns
% -------
% fz_data:
%   Inverse-transformed magnitude fMRI data [x, y, z, M].
%
% Expects
% -------
% data:
%   Denoised VST-domain data [x, y, z, M].
%
% sigma:
%   Original magnitude-domain spatial noise map [x, y, z].
%
% b:
%   Spatial inverse-VST kernel. Scalar or [bx by bz].
%
% VST_ABC:
%   Variance-stabilizing transform name or filename.
%   Default is 'B'.
%
% See also:
%   riceVST_EUI
%   perform_riceVST3
%   estimate_noise_vst3
%
% Original copyright:
% Copyright (C) 2019 CMRR at UMN
% Author: Xiaoping Wu <xpwu@cmrr.umn.edu>
%
% Low-memory modification:
% Process one patch at a time instead of allocating large 5D arrays.

%% Input preparation

if ~isa(data, 'single')
    data = single(data);
end

if ~isa(sigma, 'single')
    sigma = single(sigma);
end

[sx, sy, sz, M] = size(data);

if ~exist('VST_ABC', 'var') || isempty(VST_ABC)
    VST_ABC = 'B';
end

if nargin < 3 || isempty(b)
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

if any(b < 1) || any(mod(b, 1) ~= 0)
    error('All inverse-VST kernel dimensions must be positive integers.');
end

if bx > sx || by > sy || bz > sz
    error(['Inverse-VST kernel [%d %d %d] exceeds image dimensions ', ...
           '[%d %d %d].'], ...
           bx, by, bz, sx, sy, sz);
end

if ~isequal(size(sigma), [sx, sy, sz])
    error(['The sigma map dimensions must match the first three data ', ...
           'dimensions. Data: [%d %d %d], sigma: [%d %d %d].'], ...
           sx, sy, sz, ...
           size(sigma, 1), size(sigma, 2), size(sigma, 3));
end

% Preserve the original step length.
step = 1;

time0 = clock;

fprintf('--------start VST EUI --------\n');
fprintf('Input dimensions: %d x %d x %d x %d\n', ...
    sx, sy, sz, M);

fprintf('Inverse-VST kernel: %d x %d x %d\n', ...
    bx, by, bz);

fprintf('VST type: %s\n', VST_ABC);

%% Low-memory accumulators

% One complete 4D output accumulator.
Ys = zeros(sx, sy, sz, M, 'single');

% The overlap count is identical across time points, so a 3D map is enough.
W = zeros(sx, sy, sz, 'single');

%% Generate patch starting positions

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

disp('-> run exact unbiased inverse VST using low-memory serial processing...');

%% Process one patch at a time

for ii = 1:nPatchX

    i = indices_x(ii);

    fprintf('--- EUI x patch %d of %d, starting at x=%d ---\n', ...
        ii, nPatchX, i);

    for jj = 1:nPatchY

        j = indices_y(jj);

        for kk = 1:nPatchZ

            k = indices_z(kk);

            % Extract the current denoised VST-domain patch.
            B1 = data( ...
                i:i+bx-1, ...
                j:j+by-1, ...
                k:k+bz-1, ...
                :);

            % Extract the corresponding magnitude-noise patch.
            Sig = sigma( ...
                i:i+bx-1, ...
                j:j+by-1, ...
                k:k+bz-1);

            % Preserve the original patch-level mean sigma operation.
            patchSigma = mean(Sig(:), 'double');

            if ~isscalar(patchSigma) || ...
                    ~isfinite(patchSigma) || ...
                    patchSigma <= 0

                error(['Invalid mean sigma at inverse-VST patch start ', ...
                       '[x=%d, y=%d, z=%d].'], ...
                       i, j, k);
            end

            % Exact unbiased inverse VST.
            inversePatch = riceVST_EUI( ...
                B1, ...
                patchSigma, ...
                VST_ABC);

            inversePatch = single(inversePatch);

            if ~isequal(size(inversePatch), size(B1))
                error(['riceVST_EUI returned an unexpected output size ', ...
                       'at patch start [x=%d, y=%d, z=%d].'], ...
                       i, j, k);
            end

            if any(~isfinite(inversePatch(:)))
                error(['riceVST_EUI returned NaN or Inf at patch start ', ...
                       '[x=%d, y=%d, z=%d].'], ...
                       i, j, k);
            end

            % Accumulate the inverse-transformed patch.
            Ys( ...
                i:i+bx-1, ...
                j:j+by-1, ...
                k:k+bz-1, ...
                :) = ...
            Ys( ...
                i:i+bx-1, ...
                j:j+by-1, ...
                k:k+bz-1, ...
                :) + inversePatch;

            % Count one estimate for each spatial voxel in the patch.
            W( ...
                i:i+bx-1, ...
                j:j+by-1, ...
                k:k+bz-1) = ...
            W( ...
                i:i+bx-1, ...
                j:j+by-1, ...
                k:k+bz-1) + 1;

        end
    end
end

%% Aggregate overlapping inverse-VST estimates

disp('-> aggregate overlapping inverse-VST estimates...');

if any(W(:) == 0)
    error(['At least one voxel did not receive an inverse-VST estimate. ', ...
           'Check kernel dimensions and sliding-window indices.']);
end

fz_data = bsxfun(@rdivide, Ys, W);

if any(~isfinite(fz_data(:)))
    error('The final inverse-VST output contains NaN or Inf values.');
end

% Magnitude MRI values should not be negative.
fz_data(fz_data < 0) = 0;

fprintf('Inverse-VST output range: minimum=%g, maximum=%g\n', ...
    min(fz_data(:)), max(fz_data(:)));

fprintf('Total elapsed time = %f min\n\n', ...
    etime(clock, time0) / 60);

end
