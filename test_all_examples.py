"""
Example usage of the updated resample_nifti.py script
Demonstrates all ways to specify new spacing
"""

import numpy as np
import nibabel as nib
from resample_nifti import resample_nifti_spacing

print("="*70)
print("NIFTI RESAMPLING - USAGE EXAMPLES")
print("="*70)

# Create a dummy NIfTI file for demonstration
print("\n1. Creating a test volume (100 micron isotropic)...")
test_data = np.random.rand(100, 100, 100).astype(np.float32)
test_affine = np.array([
    [-0.1, 0.0, 0.0, 0.0],
    [0.0, -0.1, 0.0, 0.0],
    [0.0, 0.0, 0.1, 0.0],
    [0.0, 0.0, 0.0, 1.0]
])
test_img = nib.Nifti1Image(test_data, test_affine)
nib.save(test_img, 'test_100um_iso.nii.gz')
print("   Created: test_100um_iso.nii.gz (100x100x100 voxels at 100μm isotropic)")

# Example 1: Change only Z-spacing using individual parameter
print("\n" + "="*70)
print("EXAMPLE 1: Change only Z-spacing (100x100x150 μm)")
print("="*70)
resample_nifti_spacing(
    'test_100um_iso.nii.gz',
    'test_example1_z150.nii.gz',
    new_z_spacing=0.15
)

# Example 2: Change X and Z spacing
print("\n" + "="*70)
print("EXAMPLE 2: Change X and Z spacing (200x100x150 μm)")
print("="*70)
resample_nifti_spacing(
    'test_100um_iso.nii.gz',
    'test_example2_xz.nii.gz',
    new_x_spacing=0.2,
    new_z_spacing=0.15
)

# Example 3: Using tuple input
print("\n" + "="*70)
print("EXAMPLE 3: Using tuple input (150x150x200 μm)")
print("="*70)
resample_nifti_spacing(
    'test_100um_iso.nii.gz',
    'test_example3_tuple.nii.gz',
    new_spacing=(0.15, 0.15, 0.2)
)

# Example 4: Upsample in one direction (higher resolution)
print("\n" + "="*70)
print("EXAMPLE 4: Upsample Z to higher resolution (100x100x50 μm)")
print("="*70)
resample_nifti_spacing(
    'test_100um_iso.nii.gz',
    'test_example4_upsample.nii.gz',
    new_z_spacing=0.05,
    interpolation_order=3  # Cubic for upsampling
)

# Example 5: Make isotropic but coarser
print("\n" + "="*70)
print("EXAMPLE 5: Downsample to coarser isotropic (200x200x200 μm)")
print("="*70)
resample_nifti_spacing(
    'test_100um_iso.nii.gz',
    'test_example5_coarse.nii.gz',
    new_spacing=(0.2, 0.2, 0.2)
)

# Example 6: Test with rotated affine
print("\n" + "="*70)
print("EXAMPLE 6: Handling rotated affine matrix")
print("="*70)
rotated_affine = np.array([
    [-0.01045284,  0.0,        0.09945219, 0.0],
    [ 0.0,        -0.1,        0.0,        0.0],
    [ 0.09945219,  0.0,        0.01045284, 0.0],
    [ 0.0,         0.0,        0.0,        1.0]
])
rotated_img = nib.Nifti1Image(test_data, rotated_affine)
nib.save(rotated_img, 'test_rotated.nii.gz')
print("   Created rotated volume")

resample_nifti_spacing(
    'test_rotated.nii.gz',
    'test_example6_rotated.nii.gz',
    new_z_spacing=0.15
)

print("\n" + "="*70)
print("VERIFICATION")
print("="*70)

def verify_output(filename, expected_spacing_um):
    img = nib.load(filename)
    spacing = np.sqrt(np.sum(img.affine[:3, :3] ** 2, axis=0))
    spacing_um = spacing * 1000
    print(f"\n{filename}:")
    print(f"  Shape: {img.shape}")
    print(f"  Spacing (μm): x={spacing_um[0]:.1f}, y={spacing_um[1]:.1f}, z={spacing_um[2]:.1f}")
    print(f"  Expected: {expected_spacing_um}")
    
    # Check if close to expected
    expected = np.array([float(x) for x in expected_spacing_um.replace('μm', '').replace('x', ',').split(',')])
    if np.allclose(spacing_um, expected, rtol=0.01):
        print("  ✓ PASS")
    else:
        print("  ✗ FAIL")

verify_output('test_example1_z150.nii.gz', '100x100x150μm')
verify_output('test_example2_xz.nii.gz', '200x100x150μm')
verify_output('test_example3_tuple.nii.gz', '150x150x200μm')
verify_output('test_example4_upsample.nii.gz', '100x100x50μm')
verify_output('test_example5_coarse.nii.gz', '200x200x200μm')
verify_output('test_example6_rotated.nii.gz', '100x100x150μm')

print("\n" + "="*70)
print("COMMAND LINE EQUIVALENTS")
print("="*70)
print("""
The above examples can also be run from command line:

# Example 1 - Change only Z
python resample_nifti.py input.nii.gz output.nii.gz --z-spacing 0.15

# Example 2 - Change X and Z
python resample_nifti.py input.nii.gz output.nii.gz --x-spacing 0.2 --z-spacing 0.15

# Example 3 - Using tuple
python resample_nifti.py input.nii.gz output.nii.gz --spacing 0.15 0.15 0.2

# Example 4 - Upsample with cubic
python resample_nifti.py input.nii.gz output.nii.gz --z-spacing 0.05 --order 3

# Example 5 - Coarser isotropic
python resample_nifti.py input.nii.gz output.nii.gz --spacing 0.2 0.2 0.2

# Example 6 - Same as Example 1
python resample_nifti.py rotated.nii.gz output.nii.gz --z-spacing 0.15
""")

print("\n" + "="*70)
print("CLEANUP")
print("="*70)
import os
test_files = [
    'test_100um_iso.nii.gz',
    'test_rotated.nii.gz',
    'test_example1_z150.nii.gz',
    'test_example2_xz.nii.gz',
    'test_example3_tuple.nii.gz',
    'test_example4_upsample.nii.gz',
    'test_example5_coarse.nii.gz',
    'test_example6_rotated.nii.gz'
]
for f in test_files:
    if os.path.exists(f):
        os.remove(f)
        print(f"  Removed: {f}")

print("\nAll done!")
