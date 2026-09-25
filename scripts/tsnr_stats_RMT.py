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
