# NIfTI 3D Volume Resampling: Isotropic to Anisotropic or 'any'tropic

Tools for resampling MRI NIfTI volumes from 100 micron or any isotropic spacing to a new micron spacing.

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

**Basic Examples:**

```bash
# Change only z-spacing to 150 microns (keeps x,y at original spacing)
python resample_nifti.py input.nii.gz output.nii.gz --z-spacing 0.15

# Change multiple axes individually
python resample_nifti.py input.nii.gz output.nii.gz --x-spacing 0.2 --z-spacing 0.15

# Use tuple input to specify all three axes at once
python resample_nifti.py input.nii.gz output.nii.gz --spacing 0.1 0.1 0.15

# Upsample z-axis to higher resolution (50 microns) with cubic interpolation
python resample_nifti.py input.nii.gz output.nii.gz --z-spacing 0.05 --order 3

# Make volume isotropic at 200 microns
python resample_nifti.py input.nii.gz output.nii.gz --spacing 0.2 0.2 0.2

# With comparison statistics
python resample_nifti.py input.nii.gz output.nii.gz --z-spacing 0.15 --compare
```

**Advanced Examples:**

```bash
# Downsample x and y, upsample z (for highly anisotropic data)
python resample_nifti.py input.nii.gz output.nii.gz --x-spacing 0.5 --y-spacing 0.5 --z-spacing 0.05

# Convert from anisotropic to isotropic
python resample_nifti.py aniso.nii.gz iso.nii.gz --spacing 0.15 0.15 0.15 --order 1
```

### Python Script Usage

**Method 1: Individual axis parameters**

```python
from resample_nifti import resample_nifti_spacing

# Change only z-spacing (x and y keep original spacing)
resample_nifti_spacing(
    'input.nii.gz',
    'output.nii.gz',
    new_z_spacing=0.15  # 150 microns
)

# Change x and z
resample_nifti_spacing(
    'input.nii.gz',
    'output.nii.gz',
    new_x_spacing=0.2,   # 200 microns
    new_z_spacing=0.15   # 150 microns
    # y keeps original spacing
)

# Change all three axes
resample_nifti_spacing(
    'input.nii.gz',
    'output.nii.gz',
    new_x_spacing=0.15,
    new_y_spacing=0.15,
    new_z_spacing=0.2
)
```

**Method 2: Tuple input**

```python
# Specify all axes at once as (x, y, z) tuple
resample_nifti_spacing(
    'input.nii.gz',
    'output.nii.gz',
    new_spacing=(0.1, 0.1, 0.15)  # 100x100x150 microns
)
```

**Method 3: Low-level approach (for custom workflows)**

```python
import nibabel as nib
from scipy.ndimage import zoom
import numpy as np

# Load your NIfTI file
img = nib.load('your_volume.nii.gz')
data = img.get_fdata()
affine = img.affine

# Get current spacing (handles rotation properly!)
current_spacing = np.sqrt(np.sum(affine[:3, :3] ** 2, axis=0))

# Set target spacing
target_spacing = np.array([0.1, 0.1, 0.15])  # x, y, z in mm

# Calculate zoom factors (original/new)
zoom_factors = current_spacing / target_spacing

# Resample (order=1 for linear interpolation)
resampled_data = zoom(data, zoom_factors, order=1, mode='nearest')

# Update affine matrix (important for rotated volumes!)
new_affine = affine.copy()
for i in range(3):
    col_vector = affine[:3, i]
    current_magnitude = np.sqrt(np.sum(col_vector ** 2))
    if current_magnitude > 0:
        scale_factor = target_spacing[i] / current_magnitude
        new_affine[:3, i] = col_vector * scale_factor

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

### Extracting Spacing from Affine Matrix

**Important**: The voxel spacing is **not** simply the diagonal of the affine matrix when there's rotation. You must calculate the magnitude of each column vector:

```python
# WRONG method - only works for axis-aligned volumes
wrong_spacing = np.abs(affine.diagonal()[:3])

# CORRECT method - handles rotation properly
correct_spacing = np.sqrt(np.sum(affine[:3, :3] ** 2, axis=0))
```

For example, with this affine matrix:
```
[[-0.01045284  0.          0.09945219 -2.93459034]
 [ 0.         -0.1         0.          8.        ]
 [ 0.09945219  0.          0.01045284 -5.08460236]
 [ 0.          0.          0.          1.        ]]
```

- Wrong method gives: [10.5, 100.0, 10.5] μm
- Correct method gives: [100.0, 100.0, 100.0] μm ✓

### What's Happening

When resampling from 100μm isotropic to 100×100×150μm:
- **X and Y directions**: No change (zoom factor = 1.0)
- **Z direction**: Downsampling (zoom factor = 100/150 = 0.667)

This means:
- Number of slices reduces by ~33%
- File size reduces by ~33%
- Z-direction resolution decreases from 100μm to 150μm

**General Case - Multi-Axis Resampling:**

The script can resample any combination of axes:
- **Downsampling** (zoom < 1): Reduces resolution, smaller file
- **Upsampling** (zoom > 1): Increases resolution, larger file  
- **Mixed**: Different operations per axis (e.g., downsample x,y but upsample z)

Examples:
```
100μm iso → 200×100×150μm:  zoom = (0.5, 1.0, 0.667)  # x down, y same, z down
100μm iso → 50×50×100μm:    zoom = (2.0, 2.0, 1.0)    # x,y up, z same
150×150×100μm → 100μm iso:  zoom = (1.5, 1.5, 1.0)    # make isotropic
```

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
├── resample_nifti.py          # Main command-line script
├── nifti_resampling_demo.ipynb # Interactive Jupyter notebook
├── requirements.txt            # Python dependencies
└── README.md                   # This file
```

## Command Line Options

```
positional arguments:
  input                 Input NIfTI file path
  output                Output NIfTI file path

spacing arguments (mutually exclusive groups):
  --spacing X Y Z       New spacing as three values in mm (e.g., 0.1 0.1 0.15)
  --xyz-spacing X Y Z   Alternative syntax for --spacing
  
  OR use individual axes (can be combined):
  --x-spacing X         New x-spacing in mm (keeps original if not specified)
  --y-spacing Y         New y-spacing in mm (keeps original if not specified)
  --z-spacing Z         New z-spacing in mm (keeps original if not specified)

optional arguments:
  -h, --help            Show this help message and exit
  --order {0,1,2,3}     Interpolation order:
                        0=nearest, 1=linear (default), 2=quadratic, 3=cubic
  --compare             Compare original and resampled volumes
  --quiet               Suppress verbose output
```

**Important Notes:**
- You can use either `--spacing X Y Z` for all axes, OR individual `--x-spacing`, `--y-spacing`, `--z-spacing`
- If using individual parameters, unspecified axes keep their original spacing
- Spacing values are in millimeters (mm): 0.1 mm = 100 microns, 0.15 mm = 150 microns
- The script correctly handles rotated affine matrices

## Examples

### Example 1: Standard Resampling (Single Axis)
```bash
python resample_nifti.py brain_100um.nii.gz brain_150um.nii.gz --z-spacing 0.15
```

Output:
```
Loading NIfTI file: brain_100um.nii.gz

Original volume shape: (512, 512, 300)
Original spacing (mm): x=0.1000, y=0.1000, z=0.1000
Original spacing (μm): x=100.00, y=100.00, z=100.00

Changing spacing: z=0.15

Target spacing (mm): x=0.1000, y=0.1000, z=0.1500
Target spacing (μm): x=100.00, y=100.00, z=150.00

Zoom factors: x=1.0000, y=1.0000, z=0.6667
  x-axis: No change
  y-axis: No change
  z-axis: Downsampling (0.6667x)

Interpolation method: Linear (Trilinear)

Resampling volume...
Resampled volume shape: (512, 512, 200)
  X: 512 → 512 (+0.0%)
  Y: 512 → 512 (+0.0%)
  Z: 300 → 200 (-33.3%)

Verified new spacing (mm): x=0.1000, y=0.1000, z=0.1500
Verified new spacing (μm): x=100.00, y=100.00, z=150.00

Done!
```

### Example 2: Multi-Axis Resampling
```bash
python resample_nifti.py scan.nii.gz scan_aniso.nii.gz --x-spacing 0.2 --z-spacing 0.15
```

This creates a volume with 200×100×150 μm spacing.

### Example 3: Using Tuple Input
```bash
python resample_nifti.py scan.nii.gz scan_iso.nii.gz --spacing 0.15 0.15 0.15
```

Makes the volume isotropic at 150 μm.

### Example 4: Upsampling for Super-Resolution
```bash
python resample_nifti.py scan.nii.gz scan_hires.nii.gz --z-spacing 0.05 --order 3
```

Upsamples z-axis to 50 μm using cubic interpolation.

### Example 5: Batch Processing
```bash
for file in *_100um.nii.gz; do
    output="${file/_100um/_150um}"
    python resample_nifti.py "$file" "$output" --z-spacing 0.15
done
```

### Example 6: Handling Rotated Volumes
```bash
# Works correctly even with oblique acquisitions
python resample_nifti.py oblique_scan.nii.gz resampled.nii.gz --spacing 0.2 0.2 0.2
```

The script automatically handles rotation in the affine matrix.

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
