# RMT-fMRI-MEICA-pipeline

Low-memory magnitude-domain RMT denoising for multi-echo fMRI, followed by BIDS restoration, fMRIPrep echo-wise preprocessing, tedana ME-ICA, and volume-space tSNR evaluation.

This repository documents the complete workflow that was successfully used for a three-echo 7 T resting-state fMRI dataset:

```text
Magnitude multi-echo fMRI
    -> spatial Rician noise estimation
    -> forward Rice VST
    -> low-memory patch-based OP-SSVD denoising
    -> exact unbiased inverse Rice VST
    -> NIfTI header restoration
    -> BIDS dataset construction
    -> fMRIPrep echo-wise preprocessing
    -> tedana with an explicit fMRIPrep brain mask
    -> TE2, optimal-combination, and ME-ICA tSNR evaluation
```

The RMT stage is applied independently to Echo 1, Echo 2, and Echo 3 before fMRIPrep. The tested MATLAB workflow processes one echo at a time. Before each run, the input path, cache names, and output names are manually changed to the matching echo number.

## Scope

This repository focuses on volume-space processing and comparison. Surface mapping is not required for the TE2, optimal-combination, ME-ICA, tSNR, and future volume-space functional-connectivity analyses described here.

## Attribution

The original RMT, SSVD, and Rice-VST code identifies Xiaoping Wu and Wei Zhu at the Center for Magnetic Resonance Research, University of Minnesota, as original authors. This repository documents and integrates the tested low-memory serial execution, operator-norm shrinkage configuration, echo-specific caching, BIDS/header repair, fMRIPrep, tedana, and tSNR quality-control workflow. 
# 0. Repository layout

```text
RMT-fMRI-MEICA-pipeline/
|
|-- README.md
|-- requirements.txt
|
|-- matlab/
|   |-- fmriDenoising_demo.m
|   |-- fmriDenoise_SSVD2.m
|   |-- estimate_noise_vst3.m
|   |-- perform_riceVST3.m
|   |-- denoise_ssvd.m
|   |-- ssvd.m
|   |-- perform_riceVST_EUI3.m
|   |-- MCSure.m
|   `-- dependencies/
|       |-- riceVST_sigmaEst.m
|       |-- riceVST.m
|       |-- riceVST_EUI.m
|       |-- squish.m
|       |-- gunziptemp.m
|       |-- load_untouch_nii.m
|       |-- make_nii.m
|       `-- save_nii.m
|
|-- scripts/
|   |-- fix_RMT_header.py
|   |-- check_fmriprep_echoes.py
|   |-- check_mask_size.py
|   |-- compare_masks.py
|   |-- make_tsnr_fmriprepmask.py
|   |-- mask_TE2_RMT.py
|   `-- tsnr_stats_RMT.py
|
`-- examples/
    `-- expected_outputs.txt
```

Only files used in, or required by, the successful processing path should be placed in the main workflow folders. The unsuccessful MATLAB header-repair script that depended on `save_untouch_nii` is not part of the main workflow.


## Path-replacement index

Before running the workflow, replace the numbered items used throughout this README:

```text
①  Echo-specific input NIfTI path
②  RMT intermediate/output directory
③  Echo-specific sigmaVST cache filename
④  Echo-specific forward-VST cache filename
⑤  Echo-specific OP-SSVD cache filename
⑥  Echo-specific final RMT output filename
⑦  Echo-specific sigma/rank/sigmaVST output filenames
⑧  Original functional-data directory with correct NIfTI headers
⑨  RMT BIDS functional directory
⑩  RMT BIDS root mounted into fMRIPrep
⑪  fMRIPrep derivatives directory
⑫  fMRIPrep work directory
⑬  FreeSurfer license file
⑭  Participant label
⑮  Three fMRIPrep echo-wise output paths
⑯  Tedana echo-1 input
⑰  Tedana echo-2 input
⑱  Tedana echo-3 input
⑲  fMRIPrep BOLD brain mask
⑳  Final tedana output directory
㉑  Final masked-tedana directory used for tSNR
㉒  tSNR output directory
```

# 1. Environment and setup

## 1.1 Tested software

The successful downstream run used:

```text
Windows CMD
MATLAB
Python 3.12
nibabel 5.3.2
Docker 29.2.0
fMRIPrep 25.2.5
Tedana 25.1.0
```

Check the installed versions:

```cmd
python --version
python -c "import nibabel; print(nibabel.__version__)"
docker --version
tedana --version
```

## 1.2 Python packages

Example `requirements.txt`:

```text
nibabel==5.3.2
numpy
scipy
nilearn
scikit-learn
tedana==25.1.0
```

Install from Windows CMD:

```cmd
python -m pip install -r requirements.txt
```

## 1.3 Docker and fMRIPrep

Install Docker Desktop, verify Docker, and pull fMRIPrep:

```cmd
docker run --rm hello-world
docker pull nipreps/fmriprep:latest
docker images
```

## 1.4 FreeSurfer license

The tested license location was:

```text
D:\freesurfer\license.txt
```

If a license is mounted, mount the file itself:

```cmd
-v "D:/freesurfer/license.txt:/opt/freesurfer/license.txt"
```

Do not mount a local folder containing only `license.txt` over the complete container directory `/opt/freesurfer`, because that hides the FreeSurfer executables inside the container.

# 2. Input data

## 2.1 Tested three-echo input

```text
D:\MEICA\MECIA TEST 7T\sub-01\func\
    sub-01_task-rest_echo-1_bold.nii
    sub-01_task-rest_echo-1_bold.json
    sub-01_task-rest_echo-2_bold.nii
    sub-01_task-rest_echo-2_bold.json
    sub-01_task-rest_echo-3_bold.nii
    sub-01_task-rest_echo-3_bold.json
```

The tested acquisition metadata were:

```text
Echo 1: 0.01520 s
Echo 2: 0.03423 s
Echo 3: 0.05326 s
TR:     2.1 s
```

The functional NIfTI shape was:

```text
100 x 102 x 84 x 300
```

## 2.2 Important input assumption

This successful workflow is a magnitude-domain workflow. It uses spatial Rician-noise estimation, Rice variance stabilization, and exact unbiased inverse Rice VST. It is not the complex-domain, no-VST branch.

# 3. RMT denoising

## 3.1 Tested algorithm order

For each echo separately:

```text
4D magnitude NIfTI
    -> optional global rescaling if values exceed int16 range
    -> estimate_noise_vst3
    -> perform_riceVST3
    -> denoise_ssvd
    -> ssvd(..., 'op', 'ssvd', 1:10)
    -> rank-weighted overlapping-patch aggregation
    -> perform_riceVST_EUI3
    -> RMT-denoised magnitude NIfTI
```

## 3.2 Tested parameters

The final successful `fmriDenoising_demo.m` configuration used:

```matlab
ks = [7 7 3];
ws = [7 7 3];

[epis_vstssvd, rank_vstssvd, sigma_vstssvd, sigmaVST] = ...
    cellfun(@(x) fmriDenoise_SSVD2( ...
        x, ...       % 4D magnitude data
        ks, ...      % Rician noise-estimation and VST kernel
        ws, ...      % OP-SSVD patch window
        'A', ...     % forward and inverse VST type
        0, ...       % Gaussian weighting request
        1, ...       % rank-weighting request
        1, ...       % enable VST
        1:10), ...   % SSVD moment orders
        epis, ...
        'UniformOutput', 0);
```

Tested settings:

```text
Noise/VST spatial kernel:       [7 7 3]
OP-SSVD spatial patch:          [7 7 3]
Patch step:                     1
Forward/inverse VST type:       A
Noise-estimation VST type:      B in the tested fmriDenoise_SSVD2.m
SSVD rank/noise estimator:      multiple-criteria SSVD
SSVD moment orders:             1:10
Singular-value manipulation:    operator-norm optimal shrinkage
Patch aggregation weight:       1 / (1 + estimated patch rank)
VST enabled:                    yes
```

### Weighting note

In the tested OP implementation, Gaussian weighting is disabled and the whole-patch weight `1/(1+R)` is always applied. `wantGaussWeighting` and `wantRankWeighting` remain in the function signature for compatibility with older code, but they do not override the OP aggregation rule.

## 3.3 Run Echo 1, Echo 2, and Echo 3 separately

### MATLAB files used in this step

```text
matlab/
├── fmriDenoising_demo.m
├── fmriDenoise_SSVD2.m
├── estimate_noise_vst3.m
├── perform_riceVST3.m
├── denoise_ssvd.m
├── ssvd.m
├── perform_riceVST_EUI3.m
├── MCSure.m
└── dependencies/
    ├── riceVST_sigmaEst.m
    ├── riceVST.m
    ├── riceVST_EUI.m
    ├── squish.m
    ├── gunziptemp.m
    ├── load_untouch_nii.m
    ├── make_nii.m
    └── save_nii.m
```

The tested workflow processes one echo at a time.

Before running each echo, manually update:

1. The input echo path in `fmriDenoising_demo.m`.
2. The echo-specific cache filenames in `fmriDenoise_SSVD2.m`.
3. The final RMT and QC output filenames in `fmriDenoising_demo.m`.

The currently saved examples show the final Echo 3 configuration. Echo 1 and Echo 2 use the same processing parameters. Only the echo-specific input, cache, and output names are changed.

Do not reuse one echo's cache files for a different echo.

### Replace: ① to ⑦

---

### 3.3.1 File: `matlab/fmriDenoising_demo.m`

Replace the input path, intermediate directory, and output filenames.

The complete tested driver structure is shown below.

```matlab
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
```

### Echo-specific replacements in `fmriDenoising_demo.m`

#### Echo 1

```matlab
epifilenames = {
    'D:\MEICA\MECIA TEST 7T\sub-01\func\sub-01_task-rest_echo-1_bold.nii'
};
```

Use these output filenames:

```text
echo1_RMT_OP.nii
echo1_sigma.nii
echo1_rank.nii
echo1_sigmaVST.nii
```

#### Echo 2

```matlab
epifilenames = {
    'D:\MEICA\MECIA TEST 7T\sub-01\func\sub-01_task-rest_echo-2_bold.nii'
};
```

Use these output filenames:

```text
echo2_RMT_OP.nii
echo2_sigma.nii
echo2_rank.nii
echo2_sigmaVST.nii
```

#### Echo 3

```matlab
epifilenames = {
    'D:\MEICA\MECIA TEST 7T\sub-01\func\sub-01_task-rest_echo-3_bold.nii'
};
```

Use these output filenames:

```text
echo3_RMT_OP.nii
echo3_sigma.nii
echo3_rank.nii
echo3_sigmaVST.nii
```

---

### 3.3.2 File: `matlab/fmriDenoise_SSVD2.m`

This file controls:

```text
Rician noise estimation
Forward Rice VST
OP-SSVD denoising
Exact unbiased inverse Rice VST
Intermediate cache loading and saving
```

The full implementation must be placed in:

```text
matlab/fmriDenoise_SSVD2.m
```

Before each echo run, replace the three cache filenames shown below.

```matlab
%% Rician noise-estimation cache

% ③ Change this filename to match the echo being processed.
%
% Echo 1:
% echo1_sigmaVST_cache.mat
%
% Echo 2:
% echo2_sigmaVST_cache.mat
%
% Echo 3:
% echo3_sigmaVST_cache.mat

sigmaCacheFile = ...
    'D:\RMT\test result\echo3_sigmaVST_cache.mat';


%% Forward-VST cache

% ④ Change this filename to match the echo being processed.
%
% Echo 1:
% echo1_forwardVST_cache.mat
%
% Echo 2:
% echo2_forwardVST_cache.mat
%
% Echo 3:
% echo3_forwardVST_cache.mat

vstCacheFile = ...
    'D:\RMT\test result\echo3_forwardVST_cache.mat';


%% OP-SSVD cache

% ⑤ Change this filename to match the echo being processed.
%
% Echo 1:
% echo1_OPSSVD_cache.mat
%
% Echo 2:
% echo2_OPSSVD_cache.mat
%
% Echo 3:
% echo3_OPSSVD_cache.mat

# 4. Restore the NIfTI header

## 4.1 Why header restoration is required

The tested MATLAB save command used:

```matlab
save_nii(make_nii(x, episize), output_file)
```

This created a new header and did not retain the original temporal and spatial metadata. Before repair, the RMT output showed:

```text
pixdim(5) = 1
xyzt_units = 0
```

The original functional NIfTI showed:

```text
pixdim(5) = 2.1
xyzt_units = 10
```

The mismatch caused BIDS validation errors for repetition-time units and repetition-time disagreement.

## 4.2 Successful fix: `scripts/fix_RMT_header.py`

Replace: ⑧ and ⑨

```python
import nibabel as nib
import os

# ⑧ Original, unprocessed echo directory with correct headers
orig_dir = r"D:\MEICA\MECIA TEST 7T\sub-01\func"

# ⑨ RMT BIDS functional directory
rmt_dir = r"D:\RMT\RMT_BIDS\sub-01\func"

files = [
    "sub-01_task-rest_echo-1_bold.nii",
    "sub-01_task-rest_echo-2_bold.nii",
    "sub-01_task-rest_echo-3_bold.nii",
]

for fname in files:
    print("\n=================================")
    print("Processing:", fname)

    orig_file = os.path.join(orig_dir, fname)
    rmt_file = os.path.join(rmt_dir, fname)
    fixed_file = rmt_file.replace(".nii", "_FIXED.nii")

    orig_img = nib.load(orig_file)
    rmt_img = nib.load(rmt_file)

    print("Original shape:", orig_img.shape)
    print("RMT shape     :", rmt_img.shape)

    if orig_img.shape != rmt_img.shape:
        raise RuntimeError(f"Dimension mismatch: {fname}")

    fixed_img = nib.Nifti1Image(
        rmt_img.get_fdata(dtype="float32"),
        affine=orig_img.affine,
        header=orig_img.header.copy(),
    )

    nib.save(fixed_img, fixed_file)

    chk = nib.load(fixed_file)
    print("Saved:", fixed_file)
    print("TR =", chk.header["pixdim"][4])
    print("xyzt_units =", chk.header["xyzt_units"])

print("\nFinished.")
```

Run:

```cmd
python scripts\fix_RMT_header.py
```

Expected validation:

```text
TR = 2.1
xyzt_units = 10
```

After validation, move the original bad-header copies out of the BIDS directory and rename the fixed files to the BIDS-valid names:

```text
sub-01_task-rest_echo-1_bold.nii
sub-01_task-rest_echo-2_bold.nii
sub-01_task-rest_echo-3_bold.nii
```

# 5. Construct the RMT BIDS dataset

## 5.1 Tested layout

```text
D:\RMT\RMT_BIDS\
|
|-- dataset_description.json
|
`-- sub-01\
    |-- anat\
    |   |-- sub-01_acq-1_T1w.nii
    |   |-- sub-01_acq-1_T1w.json
    |   |-- sub-01_acq-2_T1w.nii
    |   |-- sub-01_acq-2_T1w.json
    |   |-- sub-01_acq-3_T1w.nii
    |   `-- sub-01_acq-3_T1w.json
    |
    |-- fmap\
    |   |-- sub-01_dir-AP_epi.nii
    |   |-- sub-01_dir-AP_epi.json
    |   |-- sub-01_dir-PA_epi.nii
    |   `-- sub-01_dir-PA_epi.json
    |
    `-- func\
        |-- sub-01_task-rest_echo-1_bold.nii
        |-- sub-01_task-rest_echo-1_bold.json
        |-- sub-01_task-rest_echo-2_bold.nii
        |-- sub-01_task-rest_echo-2_bold.json
        |-- sub-01_task-rest_echo-3_bold.nii
        `-- sub-01_task-rest_echo-3_bold.json
```

The anatomical and fieldmap files were copied from the original BIDS dataset. The original functional NIfTI files were replaced by the header-repaired RMT echoes. The original echo-specific JSON sidecars were retained because the acquisition metadata, including TE, TR, phase encoding, and slice timing, were unchanged by RMT.

## 5.2 Echo metadata sanity check

```text
Echo 1: EchoNumber 1, EchoTime 0.01520
Echo 2: EchoNumber 2, EchoTime 0.03423
Echo 3: EchoNumber 3, EchoTime 0.05326
All echoes: RepetitionTime 2.1
```

# 6. fMRIPrep echo-wise preprocessing

This stage prepares aligned echo-wise BOLD series for tedana.

## 6.1 Replace: ⑩ to ⑭

```cmd
REM ⑩ RMT BIDS root:       D:/RMT/RMT_BIDS
REM ⑪ derivatives output:  D:/RMT/fmriprep_derivatives
REM ⑫ work directory:      D:/RMT/fmriprep_work
REM ⑬ license file:        D:/freesurfer/license.txt
REM ⑭ participant label:   01

docker run --rm -it ^
  -v "D:/RMT/RMT_BIDS:/data:ro" ^
  -v "D:/RMT/fmriprep_derivatives:/out" ^
  -v "D:/RMT/fmriprep_work:/work" ^
  -v "D:/freesurfer/license.txt:/opt/freesurfer/license.txt" ^
  nipreps/fmriprep:latest ^
  /data /out participant ^
  --participant-label 01 ^
  --fs-license-file /opt/freesurfer/license.txt ^
  --work-dir /work ^
  --output-spaces T1w MNI152NLin2009cAsym ^
  --me-output-echos ^
  --nprocs 8 ^
  --omp-nthreads 8 ^
  --mem 16000
```

The tested RMT run generated and validated the required echo-wise volume outputs before the optional FreeSurfer surface stage later failed during left-hemisphere topology repair. No surface outputs were used in this repository. For the documented RMT-to-tedana workflow, the required successful products are the three echo-wise `desc-preproc_bold.nii.gz` files and the BOLD brain mask listed below.

> Important reproducibility note: the exact RMT run documented here did not complete the optional FreeSurfer surface reconstruction. The required volume-space echo outputs were already written, opened successfully, and verified to have matching dimensions. This README does not treat any FreeSurfer surface product as part of the successful pipeline.

## 6.2 Expected echo-wise outputs

```text
D:\RMT\fmriprep_derivatives\sub-01\func\
    sub-01_task-rest_echo-1_desc-preproc_bold.nii.gz
    sub-01_task-rest_echo-2_desc-preproc_bold.nii.gz
    sub-01_task-rest_echo-3_desc-preproc_bold.nii.gz
    sub-01_task-rest_desc-brain_mask.nii.gz
```

## 6.3 Sanity check: `scripts/check_fmriprep_echoes.py`

Replace: ⑮

```python
import nibabel as nib

# ⑮ Replace derivative paths
files = [
    r"D:\RMT\fmriprep_derivatives\sub-01\func\sub-01_task-rest_echo-1_desc-preproc_bold.nii.gz",
    r"D:\RMT\fmriprep_derivatives\sub-01\func\sub-01_task-rest_echo-2_desc-preproc_bold.nii.gz",
    r"D:\RMT\fmriprep_derivatives\sub-01\func\sub-01_task-rest_echo-3_desc-preproc_bold.nii.gz",
]

imgs = [nib.load(f) for f in files]
shapes = [img.shape for img in imgs]

for f, shape in zip(files, shapes):
    print(f)
    print(shape)

assert len({shape[:3] for shape in shapes}) == 1
assert len({shape[3] for shape in shapes}) == 1

print("Echo-wise outputs have matching spatial dimensions and length.")
```

Run:

```cmd
python scripts\check_fmriprep_echoes.py
```

Tested result:

```text
(100, 102, 84, 300)
(100, 102, 84, 300)
(100, 102, 84, 300)
```

# 7. Full ME-ICA with tedana and the fMRIPrep brain mask

The successful final workflow explicitly supplies the fMRIPrep brain mask to tedana. This keeps the initial analysis boundary from being estimated only from the RMT-processed first echo.

Tedana still performs its internal adaptive signal-quality screening after receiving the user-defined mask.

## 7.1 Replace: ⑯ to ⑳

```cmd
REM ⑯ echo-1 fMRIPrep output
REM ⑰ echo-2 fMRIPrep output
REM ⑱ echo-3 fMRIPrep output
REM ⑲ fMRIPrep brain mask
REM ⑳ tedana output directory

tedana -d "D:\RMT\fmriprep_derivatives\sub-01\func\sub-01_task-rest_echo-1_desc-preproc_bold.nii.gz" "D:\RMT\fmriprep_derivatives\sub-01\func\sub-01_task-rest_echo-2_desc-preproc_bold.nii.gz" "D:\RMT\fmriprep_derivatives\sub-01\func\sub-01_task-rest_echo-3_desc-preproc_bold.nii.gz" -e 0.0152 0.03423 0.05326 --mask "D:\RMT\fmriprep_derivatives\sub-01\func\sub-01_task-rest_desc-brain_mask.nii.gz" --out-dir "D:\RMT\tedana_fmriprepmask" --prefix sub-01_rest_ --convention bids --fittype curvefit --n-threads 8 --lowmem --overwrite
```

### Tested options

```text
Echo times:       0.0152 0.03423 0.05326 s
Mask:             fMRIPrep BOLD brain mask
T2*/S0 fit:       curvefit
PCA selection:    AIC default
Decision tree:    tedana_orig default
ICA:              FastICA default
Seed:             42 default
Threads:          8
Low-memory mode:  enabled
Naming:           BIDS convention
```

## 7.2 Core outputs

```text
D:\RMT\tedana_fmriprepmask\
    sub-01_rest_desc-optcom_bold.nii.gz
    sub-01_rest_desc-denoised_bold.nii.gz
    sub-01_rest_T2starmap.nii.gz
    sub-01_rest_S0map.nii.gz
    sub-01_rest_desc-adaptiveGoodSignal_mask.nii.gz
    sub-01_rest_desc-ICA_components.nii.gz
    sub-01_rest_desc-ICA_mixing.tsv
    sub-01_rest_desc-tedana_metrics.tsv
    sub-01_rest_desc-tedana_registry.json
    sub-01_rest_tedana_report.html
```

This single full tedana run already produces T2*/S0, optimal combination, PCA/ICA, component classification, and denoising. Separate `t2smap` and `ica_reclassify` runs are not part of the default successful workflow.

## 7.3 Note: successful masked workflow and the earlier voxel-dropoff diagnostic

The successful workflow above uses the fMRIPrep BOLD brain mask from the start. This run completed the full tedana workflow and produced an adaptive mask with **238262** nonzero voxels, together with the final OC and ME-ICA-denoised outputs.

Before adopting that successful masked workflow, a diagnostic tedana run was performed without an explicit `--mask`. That earlier RMT run produced visibly smaller brain coverage in the tedana adaptive mask and downstream OC/ME-ICA results.

Observed adaptive-mask nonzero voxel counts for this single test dataset:

```text
Raw adaptive mask:             246372
RMT default adaptive mask:     206166
RMT + fMRIPrep mask:           238262
```

Supplying the fMRIPrep brain mask restored 32096 voxels relative to the default RMT tedana run and recovered approximately 80% of the initial 40206-voxel coverage difference. The final RMT-plus-mask coverage was within 8110 voxels, approximately 3.3%, of the Raw adaptive-mask count.

This is a dataset-specific QC observation, not a universal expected result. Tedana continued to apply adaptive signal-quality screening. In the successful masked run, the log reported that 8849 voxels in the user-defined mask did not have good signal and were removed. Therefore, the explicit fMRIPrep mask did not force all brain-mask voxels into ME-ICA. It provided a more complete initial brain boundary before tedana's adaptive screening.

For reproducible Raw-versus-RMT comparisons, use the same mask strategy in both branches, ideally supplying each branch's corresponding fMRIPrep BOLD brain mask.

# 8. Compute tSNR

Voxelwise tSNR is calculated as:

```text
tSNR = temporal mean / sample temporal standard deviation
```

No detrending is performed in this script.

## 8.1 OC and ME-ICA tSNR: `scripts/make_tsnr_fmriprepmask.py`

Replace: ㉑ and ㉒

```python
import os
import nibabel as nib


def tsnr_4d(path_4d, out_path):
    img = nib.load(path_4d)
    dat = img.get_fdata()
    tsnr = dat.mean(axis=3) / dat.std(axis=3, ddof=1)

    os.makedirs(os.path.dirname(out_path), exist_ok=True)

    nib.Nifti1Image(
        tsnr,
        img.affine,
        img.header,
    ).to_filename(out_path)

    print("Saved:", out_path)


# ㉑ Final tedana directory generated with the fMRIPrep mask
td_dir = r"D:\RMT\tedana_fmriprepmask"

# ㉒ tSNR output directory
outdir = r"D:\RMT\tSNR"
os.makedirs(outdir, exist_ok=True)

tsnr_4d(
    os.path.join(td_dir, "sub-01_rest_desc-optcom_bold.nii.gz"),
    os.path.join(outdir, "tsnr_OC_RMT_MASKED.nii.gz"),
)

tsnr_4d(
    os.path.join(td_dir, "sub-01_rest_desc-denoised_bold.nii.gz"),
    os.path.join(outdir, "tsnr_MEICA_RMT_MASKED.nii.gz"),
)

print("\nFinished.")
```

Run:

```cmd
python scripts\make_tsnr_fmriprepmask.py
```

## 8.2 TE2 tSNR

TE2 is not regenerated by tedana. It remains the fMRIPrep echo-2 output:

```text
D:\RMT\fmriprep_derivatives\sub-01\func\
sub-01_task-rest_echo-2_desc-preproc_bold.nii.gz
```

Calculate TE2 tSNR with the same `tsnr_4d` function and save:

```text
D:\RMT\tSNR\tsnr_TE2_RMT.nii.gz
```

## 8.3 Remove the TE2 brain-exterior plane: `scripts/mask_TE2_RMT.py`

```python
import nibabel as nib

# Replace these paths if needed.
tsnr_file = r"D:\RMT\tSNR\tsnr_TE2_RMT.nii.gz"
mask_file = (
    r"D:\RMT\fmriprep_derivatives\sub-01\func"
    r"\sub-01_task-rest_desc-brain_mask.nii.gz"
)
out_file = r"D:\RMT\tSNR\tsnr_TE2_RMT_MASKED.nii.gz"

tsnr_img = nib.load(tsnr_file)
mask_img = nib.load(mask_file)

tsnr = tsnr_img.get_fdata()
mask = mask_img.get_fdata()

tsnr_masked = tsnr * (mask > 0)

nib.save(
    nib.Nifti1Image(
        tsnr_masked,
        tsnr_img.affine,
        tsnr_img.header,
    ),
    out_file,
)

print("Saved:", out_file)
```

Run:

```cmd
python scripts\mask_TE2_RMT.py
```

# 9. Calculate brain-mask tSNR summary statistics

## 9.1 `scripts/tsnr_stats_RMT.py`

```python
import numpy as np
import nibabel as nib

mask_path = (
    r"D:\RMT\fmriprep_derivatives\sub-01\func"
    r"\sub-01_task-rest_desc-brain_mask.nii.gz"
)

tsnr_paths = {
    "TE2": r"D:\RMT\tSNR\tsnr_TE2_RMT.nii.gz",
    "OC": r"D:\RMT\tSNR\tsnr_OC_RMT_MASKED.nii.gz",
    "MEICA": r"D:\RMT\tSNR\tsnr_MEICA_RMT_MASKED.nii.gz",
}

mask = nib.load(mask_path).get_fdata() > 0

print("Average tSNR within brain mask:\n")

for name, path in tsnr_paths.items():
    tsnr = nib.load(path).get_fdata()
    vals = tsnr[mask]
    vals = vals[np.isfinite(vals)]
    vals = vals[vals > 0]

    print(
        f"{name}: "
        f"mean = {vals.mean():.2f}, "
        f"median = {np.median(vals):.2f}, "
        f"95% = {np.percentile(vals, 95):.2f}, "
        f"max = {vals.max():.2f}, "
        f"nVox = {vals.size}"
    )
```

Run:

```cmd
python scripts\tsnr_stats_RMT.py
```

## 9.2 Tested RMT results

```text
TE2:
mean = 17.93
median = 15.78
95th percentile = 40.61
max = 133.54
nVox = 247109

OC:
mean = 24.73
median = 22.34
95th percentile = 51.40
max = 148.43
nVox = 238262

MEICA:
mean = 69.91
median = 66.83
95th percentile = 131.16
max = 347.41
nVox = 238262
```

These values are dataset-specific examples and should not be treated as expected values for other acquisitions.

# 10. Volume-space visualization

Load the following maps in Connectome Workbench, FSLeyes, MRIcroGL, or another NIfTI viewer:

```text
D:\RMT\tSNR\tsnr_TE2_RMT_MASKED.nii.gz
D:\RMT\tSNR\tsnr_OC_RMT_MASKED.nii.gz
D:\RMT\tSNR\tsnr_MEICA_RMT_MASKED.nii.gz
```

Use the same palette and the same fixed numerical range across Raw and RMT branches. Do not compare images using independent automatic color ranges.

# 11. Successful workflow checklist

```text
[ ] Echo 1 input, cache, and output names were changed to echo1
[ ] Echo 2 input, cache, and output names were changed to echo2
[ ] Echo 3 input, cache, and output names were changed to echo3
[ ] All three RMT NIfTI files have shape 100 x 102 x 84 x 300
[ ] Original affine/header restored to each RMT echo
[ ] TR is 2.1 s after header repair
[ ] xyzt_units is 10 after header repair
[ ] BIDS validator has no TR-unit or TR-mismatch errors
[ ] fMRIPrep wrote all three desc-preproc echo files
[ ] All fMRIPrep echoes have identical shape and length
[ ] Tedana used the fMRIPrep brain mask
[ ] Tedana completed and wrote OC plus denoised BOLD
[ ] Adaptive-mask coverage was visually inspected
[ ] TE2 tSNR was masked for visualization
[ ] Raw and RMT statistics use consistent mask logic
```

# 12. Troubleshooting

## 12.1 BIDS repetition-time errors after RMT

Symptoms:

```text
REPETITION_TIME_UNITS
REPETITION_TIME_MISMATCH
```

Cause:

```text
make_nii created a new header with pixdim(5)=1 and xyzt_units=0
```

Fix:

```text
Run fix_RMT_header.py and confirm TR=2.1 and xyzt_units=10.
Do not permanently bypass this issue with --skip-bids-validation.
```

## 12.2 `mri_convert` not found in the fMRIPrep container

Cause:

```text
A local folder containing only license.txt was mounted over /opt/freesurfer.
```

Fix:

```cmd
-v "D:/freesurfer/license.txt:/opt/freesurfer/license.txt"
```

## 12.3 Tedana PCA memory error

Observed error:

```text
Unable to allocate 430 MiB
```

Successful fix:

```text
--lowmem
```

The final run retained:

```text
--n-threads 8
```

and completed the complete PCA, ICA, classification, denoising, and report workflow.

## 12.4 RMT tedana result has visibly reduced brain coverage

Use an explicit fMRIPrep BOLD brain mask:

```text
--mask sub-01_task-rest_desc-brain_mask.nii.gz
```

Then compare the final adaptive mask and denoised output visually and quantitatively. In the tested dataset, the adaptive mask increased from 206166 to 238262 nonzero voxels.

# 13. References and software documentation

- [fMRIPrep documentation](https://fmriprep.org/en/25.2.5/)
- [tedana 25.1.0 command-line documentation](https://tedana.readthedocs.io/en/25.1.0/usage.html)
- [BIDS specification](https://bids-specification.readthedocs.io/)
- [NiBabel documentation](https://nipy.org/nibabel/)

Add the original publications associated with the supplied RMT, SSVD, operator-norm shrinkage, Rice-VST, and exact unbiased inverse-transform implementations before publication or redistribution.

# 14. Reproducibility summary

The tested comparison keeps the downstream fMRIPrep and tedana settings fixed while inserting magnitude-domain VST-OP-SSVD before fMRIPrep:

```text
Raw branch:
Raw echoes -> fMRIPrep -> tedana -> TE2 / OC / ME-ICA

RMT branch:
Raw echoes -> VST-OP-SSVD -> header restoration -> fMRIPrep -> tedana -> TE2 / OC / ME-ICA
```

For formal Raw-versus-RMT comparisons, use matching fMRIPrep versions, tedana versions, output spaces, brain-mask strategy, T2* fitting method, PCA criterion, decision tree, ICA seed, tSNR definition, and downstream functional-connectivity procedure.
