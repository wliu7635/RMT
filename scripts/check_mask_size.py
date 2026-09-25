import nibabel as nib

m = nib.load(
r"D:\RMT\tedana_fmriprepmask\sub-01_rest_desc-adaptiveGoodSignal_mask.nii.gz"
).get_fdata()

print(
    "voxels =",
    (m > 0).sum()
)
