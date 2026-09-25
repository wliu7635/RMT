import nibabel as nib

# ㉓ Raw tedana adaptive mask

old_mask = nib.load(
r"D:\MEICA\MECIA-derivatives_nofs\tedana\sub-01_task-rest\sub-01_task-rest_desc-adaptiveGoodSignal_mask.nii.gz"
).get_fdata()

# ㉔ Initial RMT tedana adaptive mask

new_mask = nib.load(
r"D:\RMT\tedana\sub-01_rest_desc-adaptiveGoodSignal_mask.nii.gz"
).get_fdata()

print(
    "OLD mask voxels =",
    (old_mask > 0).sum()
)

print(
    "NEW mask voxels =",
    (new_mask > 0).sum()
)
