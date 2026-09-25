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
