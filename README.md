# NIfTI Volume Resampling: Isotropic to Anisotropic

Tools for resampling MRI NIfTI volumes from 100 micron isotropic spacing to 100×100×150 micron anisotropic spacing.

## Installation

```bash
pip install -r requirements.txt --break-system-packages
```

Or install packages individually:

```bash
pip install numpy nibabel scipy matplotlib SimpleITK --break-system-packages
```

## Quick Start

### Command Line Usage

```bash
# Basic usage with default settings (linear interpolation, 150 micron z-spacing)
python resample_nifti.py input.nii.gz output.nii.gz

# With custom z-spacing (200 microns = 0.2 mm)
python resample_nifti.py input.nii.gz output.nii.gz --z-spacing 0.20

# Using cubic interpolation
python resample_nifti.py input.nii.gz output.nii.gz --order 3

# With comparison statistics
python resample_nifti.py input.nii.gz output.nii.gz --compare
```

### Python Script Usage

```python
import nibabel as nib
from scipy.ndimage import zoom
import numpy as np

# Load your NIfTI file
img = nib.load('your_volume.nii.gz')
data = img.get_fdata()
affine = img.affine

# Get current spacing
current_spacing = np.abs(affine.diagonal()[:3])

# Calculate zoom factor (original/new)
new_z_spacing = 0.15  # 150 microns in mm
zoom_factors = [1.0, 1.0, current_spacing[2] / new_z_spacing]

# Resample (order=1 for linear interpolation)
resampled_data = zoom(data, zoom_factors, order=1, mode='nearest')

# Update affine matrix
new_affine = affine.copy()
new_affine[2, 2] = new_z_spacing if affine[2, 2] > 0 else -new_z_spacing

# Save
resampled_img = nib.Nifti1Image(resampled_data, new_affine)
nib.save(resampled_img, 'output.nii.gz')
```

### Jupyter Notebook

Open `nifti_resampling_demo.ipynb` for an interactive demonstration with visualizations.

## Interpolation Methods

The interpolation order parameter controls how new voxel values are calculated:

- **0 = Nearest Neighbor**: Fastest, but blocky. Good for label maps.
- **1 = Linear (Trilinear)**: **Recommended**. Good balance of quality and speed.
- **2 = Quadratic**: Smoother than linear, but slower.
- **3 = Cubic**: Smoothest, but can introduce ringing artifacts.

## Understanding the Transform

### What's Happening

When resampling from 100μm isotropic to 100×100×150μm:
- **X and Y directions**: No change (zoom factor = 1.0)
- **Z direction**: Downsampling (zoom factor = 100/150 = 0.667)

This means:
- Number of slices reduces by ~33%
- File size reduces by ~33%
- Z-direction resolution decreases from 100μm to 150μm

### Effects on Data Quality

**Information Loss**:
- You're going from higher to lower resolution in Z
- Fine structures < 150μm in Z may be lost or blurred
- This is **irreversible** - upsampling later won't recover lost detail

**Interpolation Artifacts**:
- **Linear**: Slight smoothing, minimal artifacts (recommended)
- **Nearest**: Blocky, stair-step artifacts in Z
- **Cubic**: Potential ringing around sharp edges

**Anisotropic Voxels**:
- Some analysis algorithms assume isotropic voxels
- Edge detection and gradient calculations may need adjustment
- Consider this when using segmentation or registration tools

## File Structure

```
.
├── resample_nifti.py           # Main command-line script
├── nifti_resampling_demo.ipynb # Interactive Jupyter notebook
├── requirements.txt            # Python dependencies
└── README.md                   # This file
```

## Command Line Options

```
positional arguments:
  input                 Input NIfTI file path
  output                Output NIfTI file path

optional arguments:
  -h, --help            show this help message and exit
  --z-spacing Z_SPACING
                        New z-spacing in mm (default: 0.15 = 150 microns)
  --order {0,1,2,3}     Interpolation order:
                        0=nearest, 1=linear (default), 2=quadratic, 3=cubic
  --compare             Compare original and resampled volumes
  --quiet               Suppress verbose output
```

## Examples

### Example 1: Standard Resampling
```bash
python resample_nifti.py brain_100um.nii.gz brain_150um.nii.gz
```

Output:
```
Loading NIfTI file: brain_100um.nii.gz

Original volume shape: (512, 512, 300)
Original spacing (mm): x=0.100, y=0.100, z=0.100
Original spacing (μm): x=100.0, y=100.0, z=100.0

Target z-spacing: 150.0 μm (0.150 mm)
Zoom factors: x=1.000, y=1.000, z=0.667
Interpolation method: Linear (Trilinear)

Resampling volume...
Resampled volume shape: (512, 512, 200)
Size reduction in z: 300 → 200 slices (33.3% reduction)

New spacing (mm): x=0.100, y=0.100, z=0.150
New spacing (μm): x=100.0, y=100.0, z=150.0

Saving resampled volume to: brain_150um.nii.gz
Done!
```

### Example 2: Different Z-spacing with Cubic Interpolation
```bash
python resample_nifti.py scan.nii.gz scan_200um.nii.gz --z-spacing 0.20 --order 3
```

### Example 3: Batch Processing
```bash
for file in *_100um.nii.gz; do
    output="${file/_100um/_150um}"
    python resample_nifti.py "$file" "$output"
done
```

## Using SimpleITK (Alternative)

For more robust handling of medical images with complex orientations:

```python
import SimpleITK as sitk

image = sitk.ReadImage('input.nii.gz')
resampler = sitk.ResampleImageFilter()
resampler.SetOutputSpacing([0.1, 0.1, 0.15])  # mm
resampler.SetSize([
    image.GetSize()[0],
    image.GetSize()[1],
    int(image.GetSize()[2] * 0.667)
])
resampler.SetInterpolator(sitk.sitkLinear)
resampler.SetOutputDirection(image.GetDirection())
resampler.SetOutputOrigin(image.GetOrigin())

resampled = resampler.Execute(image)
sitk.WriteImage(resampled, 'output.nii.gz')
```

## Quality Checks

After resampling, verify your output:

1. **Check dimensions**: Confirm Z slices reduced by expected amount
2. **Check spacing**: Use `nibabel` or `fslinfo` to verify voxel dimensions
3. **Visual inspection**: Look for unexpected artifacts
4. **Intensity statistics**: Compare min/max/mean/std with original

```python
import nibabel as nib

img = nib.load('output.nii.gz')
print(f"Shape: {img.shape}")
print(f"Spacing: {img.header['pixdim'][1:4]} mm")
print(f"Data range: [{img.get_fdata().min()}, {img.get_fdata().max()}]")
```

## Troubleshooting

**Memory errors**: Large volumes may require chunked processing or more RAM

**Unexpected artifacts**: Try different interpolation orders; cubic can cause ringing

**Wrong orientation**: SimpleITK handles orientations better than scipy for complex cases

**File won't load**: Ensure file isn't corrupted; check with `nibabel.load()` directly

## References

- [NiBabel Documentation](https://nipy.org/nibabel/)
- [SciPy ndimage Documentation](https://docs.scipy.org/doc/scipy/reference/ndimage.html)
- [SimpleITK Documentation](https://simpleitk.readthedocs.io/)

## License

This code is provided as-is for educational and research purposes.
