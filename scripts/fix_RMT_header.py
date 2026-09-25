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
