function fz_data = perform_riceVST3(data, sigma, b, VST_ABC)
% PERFORM_RICEVST3
% Apply the Rice variance-stabilizing transformation for spatially varying
% noise levels using overlapping spatial patches.
%
% This is a low-memory serial implementation of the original
% perform_riceVST3.m.
%
% The original algorithm is preserved:
%   1. Use the same overlapping sliding-window locations.
%   2. Calculate the mean noise level within each spatial patch.
%   3. Apply riceVST using that patch-level mean noise estimate.
%   4. Aggregate multiple transformed estimates for each voxel.
%   5. Divide by the number of contributing patches.
%
% The memory-intensive five-dimensional arrays data0, sigma0, Ys0, and W0
% are not created.
%
% Usage:
%   fz_data = perform_riceVST3(data, sigma, b, VST_ABC)
%
% Returns
% -------
% fz_data:
%   Variance-stabilized image series with dimensions [x, y, z, M].
%
% Expects
% -------
% data:
%   Magnitude MR image series with dimensions [x, y, z, M].
%
% sigma:
%   Spatial noise map with dimensions [x, y, z].
%
% b:
%   Spatial VST kernel. A scalar creates a b-by-b-by-b kernel.
%   A three-element vector specifies [bx, by, bz].
%
% VST_ABC:
%   Variance-stabilizing transform name or filename.
%   Default is 'B'.
%
% See also:
%   estimate_noise_vst3
%   riceVST
%   perform_riceVST_EUI3
%   riceVST_EUI
%
% Original copyright:
% Copyright (C) 2019 CMRR at UMN
% Author: Xiaoping Wu <xpwu@cmrr.umn.edu>
%
% Low-memory execution modification:
% Processes one patch at a time and stores only one 4D output accumulator
% plus one 3D overlap-count map.

%% Input preparation

% Use single precision to keep the full 4D accumulator within memory.
% The original function also converted integer input to single precision.
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

% Convert scalar kernel size to [bx by bz].
if isscalar(b)
    b = repmat(b, 1, 3);
end

% Validate kernel specification.
if numel(b) ~= 3
    error('b must be a scalar or a three-element vector [bx by bz].');
end

b = double(b(:).');

bx = b(1);
by = b(2);
bz = b(3);

if any(b < 1) || any(mod(b, 1) ~= 0)
    error('All kernel dimensions must be positive integers.');
end

if bx > sx || by > sy || bz > sz
    error(['The requested VST kernel [%d %d %d] exceeds the input ', ...
           'image dimensions [%d %d %d].'], ...
           bx, by, bz, sx, sy, sz);
end

if ~isequal(size(sigma), [sx, sy, sz])
    error(['The sigma map dimensions must match the first three data ', ...
           'dimensions. Data: [%d %d %d], sigma: [%d %d %d].'], ...
           sx, sy, sz, ...
           size(sigma, 1), size(sigma, 2), size(sigma, 3));
end

% Preserve the original sliding-window step.
step = 1;

time0 = clock;

fprintf('--------start VST --------\n');
fprintf('Input dimensions: %d x %d x %d x %d\n', sx, sy, sz, M);
fprintf('VST kernel: %d x %d x %d\n', bx, by, bz);
fprintf('VST type: %s\n', VST_ABC);

%% Allocate low-memory accumulators

% One 4D accumulator for transformed data.
Ys = zeros(sx, sy, sz, M, 'single');

% The original code allocated a 4D W array, although the patch count is
% identical for all time points. A 3D count map is sufficient.
W = zeros(sx, sy, sz, 'single');

%% Generate patch starting positions

indices_x = [1:step:sx-bx, sx-bx+1];
indices_y = [1:step:sy-by, sy-by+1];
indices_z = [1:step:sz-bz, sz-bz+1];

% Avoid duplicate final positions.
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

disp('-> run VST using low-memory serial processing...');

%% Process one spatial patch at a time

for ii = 1:nPatchX

    i = indices_x(ii);

    fprintf('--- VST x patch %d of %d, starting at x=%d ---\n', ...
        ii, nPatchX, i);

    for jj = 1:nPatchY

        j = indices_y(jj);

        for kk = 1:nPatchZ

            k = indices_z(kk);

            % Extract the current 4D patch.
            B1 = data( ...
                i:i+bx-1, ...
                j:j+by-1, ...
                k:k+bz-1, ...
                :);

            % Extract the corresponding 3D noise-map patch.
            Sig = sigma( ...
                i:i+bx-1, ...
                j:j+by-1, ...
                k:k+bz-1);

            % Use the average noise level within the current patch.
            patchSigma = mean(Sig(:), 'double');

            if ~isscalar(patchSigma) || ...
                    ~isfinite(patchSigma) || ...
                    patchSigma <= 0

                error(['Invalid mean sigma at patch start ', ...
                       '[x=%d, y=%d, z=%d].'], ...
                       i, j, k);
            end

            % Apply the same Rice VST operation used by the original code.
            rimavst = riceVST(B1, patchSigma, VST_ABC);

            % Store the transformed patch in single precision to control
            % memory usage.
            rimavst = single(rimavst);

            if ~isequal(size(rimavst), size(B1))
                error(['riceVST returned an unexpected output size at ', ...
                       'patch start [x=%d, y=%d, z=%d].'], ...
                       i, j, k);
            end

            if any(~isfinite(rimavst(:)))
                error(['riceVST returned NaN or Inf at patch start ', ...
                       '[x=%d, y=%d, z=%d].'], ...
                       i, j, k);
            end

            % Accumulate this transformed patch.
            Ys( ...
                i:i+bx-1, ...
                j:j+by-1, ...
                k:k+bz-1, ...
                :) = ...
            Ys( ...
                i:i+bx-1, ...
                j:j+by-1, ...
                k:k+bz-1, ...
                :) + rimavst;

            % Record one contribution per covered spatial voxel.
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

%% Aggregate overlapping transformed estimates

disp('-> aggregate overlapping VST estimates...');

if any(W(:) == 0)
    error(['At least one voxel did not receive a VST estimate. ', ...
           'Check the kernel dimensions and sliding-window indices.']);
end

% Divide each time point by the corresponding 3D overlap count.
% BSXFUN is used for compatibility with MATLAB versions that do not support
% implicit expansion for a 4D array divided by a 3D array.
fz_data = bsxfun(@rdivide, Ys, W);

if any(~isfinite(fz_data(:)))
    error('The final VST output contains NaN or Inf values.');
end

fprintf('VST output range: minimum=%g, maximum=%g\n', ...
    min(fz_data(:)), max(fz_data(:)));

fprintf('Total elapsed time = %f min\n\n', ...
    etime(clock, time0) / 60);

end