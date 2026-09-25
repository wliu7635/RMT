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
