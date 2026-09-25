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
