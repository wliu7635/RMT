clear;
close all;
clc;

% This driver runs magnitude-domain RMT denoising for one echo.
%
% Run this file separately for Echo 1, Echo 2, and Echo 3.
%
% Before each run, update:
%   ① input echo NIfTI
%   ② intermediate/output directory
%   ⑥ final RMT output filename
%   ⑦ sigma, rank, and sigmaVST output filenames
%
% The cache filenames inside fmriDenoise_SSVD2.m must also be
% changed to the same echo number before each run.

%% Load EPI data

% ① Input echo NIfTI
%
% Echo 1:
% epifilenames = {
%     'D:\MEICA\MECIA TEST 7T\sub-01\func\sub-01_task-rest_echo-1_bold.nii'
% };
%
% Echo 2:
% epifilenames = {
%     'D:\MEICA\MECIA TEST 7T\sub-01\func\sub-01_task-rest_echo-2_bold.nii'
% };
%
% Echo 3:
epifilenames = {
    'D:\MEICA\MECIA TEST 7T\sub-01\func\sub-01_task-rest_echo-3_bold.nii'
};

nData = length(epifilenames);

% ② Intermediate and output directory
epiIntermPath = repmat( ...
    {'D:\RMT\test result'}, ...
    1, ...
    nData);

if ~exist('D:\RMT\test result', 'dir')
    mkdir('D:\RMT\test result');
end

fprintf('Loading functional data...\n');

epis = cell(1, nData);

for p = 1:nData

    nii = load_untouch_nii( ...
        gunziptemp(epifilenames{p}));

    % Apply one global scale factor only if the data exceed the
    % signed 16-bit integer range.
    if min(nii.img(:)) < intmin('int16') || ...
            max(nii.img(:)) > intmax('int16')

        scaleFactor = ceil( ...
            max(nii.img(:)) / ...
            double(intmax('int16')));

        epis{p} = nii.img ./ scaleFactor;

        fprintf( ...
            'Applied global scale factor: %g\n', ...
            scaleFactor);
    else
        epis{p} = nii.img;
        scaleFactor = 1;

        fprintf( ...
            'No intensity scaling was required.\n');
    end

    if p == 1
        episize = nii.hdr.dime.pixdim(2:4);
        epidim = nii.hdr.dime.dim(2:5);
        epitr = 2.1;
    end
end

fprintf('Done loading EPI data.\n');

%% RMT denoising parameters

fprintf('Denoising raw data...\n');

% Noise-estimation and VST spatial kernel
ks = [7 7 3];

% OP-SSVD spatial patch window
ws = [7 7 3];

% Tested RMT configuration:
%
% VST_ABC             = 'A'
% wantGaussWeighting  = 0
% wantRankWeighting   = 1
% wantVST             = 1
% k0                  = 1:10
%
% The current OP implementation always uses whole-patch
% 1/(1 + estimated rank) weighting.

[epis_vstssvd, ...
 rank_vstssvd, ...
 sigma_vstssvd, ...
 sigmaVST] = ...
    cellfun( ...
        @(x) fmriDenoise_SSVD2( ...
            x, ...
            ks, ...
            ws, ...
            'A', ...
            0, ...
            1, ...
            1, ...
            1:10), ...
        epis, ...
        'UniformOutput', ...
        0);

%% Save the RMT outputs

% IMPORTANT:
%
% Change all filenames below to match the echo being processed.
%
% Echo 1:
%   echo1_RMT_OP.nii
%   echo1_sigma.nii
%   echo1_rank.nii
%   echo1_sigmaVST.nii
%
% Echo 2:
%   echo2_RMT_OP.nii
%   echo2_sigma.nii
%   echo2_rank.nii
%   echo2_sigmaVST.nii
%
% Echo 3:
%   echo3_RMT_OP.nii
%   echo3_sigma.nii
%   echo3_rank.nii
%   echo3_sigmaVST.nii

% ⑥ Final RMT-denoised 4D output
cellfun( ...
    @(x, y) save_nii( ...
        make_nii(x, episize), ...
        fullfile(y, 'echo3_RMT_OP.nii')), ...
    epis_vstssvd, ...
    epiIntermPath, ...
    'UniformOutput', ...
    0);

% ⑦ OP-SSVD noise-estimate output
cellfun( ...
    @(x, y) save_nii( ...
        make_nii(x, episize), ...
        fullfile(y, 'echo3_sigma.nii')), ...
    sigma_vstssvd, ...
    epiIntermPath, ...
    'UniformOutput', ...
    0);

% ⑦ OP-SSVD patch-rank output
cellfun( ...
    @(x, y) save_nii( ...
        make_nii(x, episize), ...
        fullfile(y, 'echo3_rank.nii')), ...
    rank_vstssvd, ...
    epiIntermPath, ...
    'UniformOutput', ...
    0);

% ⑦ Spatial Rician-noise map used by the VST
cellfun( ...
    @(x, y) save_nii( ...
        make_nii(x, episize), ...
        fullfile(y, 'echo3_sigmaVST.nii')), ...
    sigmaVST, ...
    epiIntermPath, ...
    'UniformOutput', ...
    0);

fprintf('RMT denoising and output writing completed.\n');
