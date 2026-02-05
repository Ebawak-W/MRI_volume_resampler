# NIfTI Resampling - Quick Reference Guide

## Basic Syntax

```bash
python resample_nifti.py INPUT.nii.gz OUTPUT.nii.gz [OPTIONS]
```

## Common Use Cases

### 1. Change Z-spacing Only (Most Common)
```bash
# 100μm iso → 100×100×150μm
python resample_nifti.py input.nii.gz output.nii.gz --z-spacing 0.15
```

### 2. Make Isotropic
```bash
# Any spacing → 150×150×150μm isotropic
python resample_nifti.py input.nii.gz output.nii.gz --spacing 0.15 0.15 0.15
```

### 3. Change Multiple Axes
```bash
# 100μm iso → 200×100×150μm
python resample_nifti.py input.nii.gz output.nii.gz --x-spacing 0.2 --z-spacing 0.15
```

### 4. Upsample for Higher Resolution
```bash
# 100μm iso → 100×100×50μm (2x resolution in Z)
python resample_nifti.py input.nii.gz output.nii.gz --z-spacing 0.05 --order 3
```

### 5. Downsample for Smaller Files
```bash
# 100μm iso → 200×200×200μm (1/8 the voxels!)
python resample_nifti.py input.nii.gz output.nii.gz --spacing 0.2 0.2 0.2
```

## Interpolation Methods

| Order | Method              | When to Use                          |
|-------|---------------------|--------------------------------------|
| 0     | Nearest Neighbor    | Label/segmentation maps only         |
| 1     | Linear (default)    | Most cases - good balance            |
| 2     | Quadratic           | Rarely needed                        |
| 3     | Cubic               | Upsampling, when quality is critical |

**Rules of thumb:**
- **Downsampling**: Use linear (order 1)
- **Upsampling**: Use cubic (order 3) for better quality
- **Label maps**: Always use nearest (order 0)

## Spacing Conversion Table

| Microns (μm) | Millimeters (mm) | Command Argument |
|--------------|------------------|------------------|
| 50           | 0.05             | `0.05`           |
| 100          | 0.1              | `0.1`            |
| 150          | 0.15             | `0.15`           |
| 200          | 0.2              | `0.2`            |
| 250          | 0.25             | `0.25`           |
| 500          | 0.5              | `0.5`            |
| 1000 (1mm)   | 1.0              | `1.0`            |

## Two Ways to Specify Spacing

### Method 1: Individual Axes (Flexible)
```bash
# Only change what you specify
python resample_nifti.py in.nii.gz out.nii.gz --z-spacing 0.15
python resample_nifti.py in.nii.gz out.nii.gz --x-spacing 0.2 --y-spacing 0.2
python resample_nifti.py in.nii.gz out.nii.gz --x-spacing 0.2 --z-spacing 0.15
```

### Method 2: Tuple (All Three at Once)
```bash
# Must specify all three: X Y Z
python resample_nifti.py in.nii.gz out.nii.gz --spacing 0.1 0.1 0.15
python resample_nifti.py in.nii.gz out.nii.gz --spacing 0.2 0.2 0.2
```

## Python Usage

### Quick Example
```python
from resample_nifti import resample_nifti_spacing

# Change only Z
resample_nifti_spacing('in.nii.gz', 'out.nii.gz', new_z_spacing=0.15)

# Use tuple
resample_nifti_spacing('in.nii.gz', 'out.nii.gz', new_spacing=(0.1, 0.1, 0.15))

# Change multiple
resample_nifti_spacing('in.nii.gz', 'out.nii.gz', 
                       new_x_spacing=0.2, new_z_spacing=0.15)
```

## File Size Impact

Resampling changes the number of voxels, which affects file size:

```
Original: 100×100×100 voxels at 100μm = 1,000,000 voxels

→ 100×100×150μm:  1,000,000 × (100/100) × (100/100) × (100/150) = 666,667 voxels (-33%)
→ 200×200×200μm:  1,000,000 × (100/200) × (100/200) × (100/200) = 125,000 voxels (-87%)
→ 50×50×50μm:     1,000,000 × (100/50)  × (100/50)  × (100/50)  = 8,000,000 voxels (+700%)
```

## Common Workflows

### Workflow 1: Reduce File Size for Faster Processing
```bash
# Downsample to half resolution (1/8 the voxels)
python resample_nifti.py large.nii.gz small.nii.gz --spacing 0.2 0.2 0.2
```

### Workflow 2: Match Spacing to Another Volume
```bash
# If target volume is 0.3×0.3×0.5mm
python resample_nifti.py source.nii.gz matched.nii.gz --spacing 0.3 0.3 0.5
```

### Workflow 3: Prepare for Algorithm Requiring Isotropic Data
```bash
# Many segmentation tools need isotropic voxels
python resample_nifti.py aniso.nii.gz iso.nii.gz --spacing 0.15 0.15 0.15 --order 1
```

### Workflow 4: Increase Z-resolution Using Super-Resolution Model
```bash
# First step: interpolate to target spacing
python resample_nifti.py low_res_z.nii.gz interp.nii.gz --z-spacing 0.05 --order 3
# Then: apply your super-resolution model to interp.nii.gz
```

## Verification

After resampling, always verify:

```python
import nibabel as nib
import numpy as np

img = nib.load('output.nii.gz')
spacing = np.sqrt(np.sum(img.affine[:3, :3] ** 2, axis=0))
print(f"Shape: {img.shape}")
print(f"Spacing (mm): {spacing}")
print(f"Spacing (μm): {spacing * 1000}")
```

Or use command line:
```bash
python resample_nifti.py input.nii.gz output.nii.gz --z-spacing 0.15 --compare
```

## Troubleshooting

**Problem**: Wrong output spacing
- **Cause**: Input had rotated affine matrix
- **Solution**: Script handles this automatically now (post-fix)

**Problem**: File is too large after upsampling
- **Cause**: Creating many more voxels than original
- **Solution**: This is expected; consider if you really need that resolution

**Problem**: Image looks blurry
- **Cause**: Downsampling with linear interpolation smooths data
- **Solution**: This is expected; you're losing information when downsampling

**Problem**: "Stair-step" artifacts
- **Cause**: Using nearest neighbor interpolation on continuous data
- **Solution**: Use `--order 1` (linear) instead

**Problem**: Ringing artifacts around edges
- **Cause**: Using cubic interpolation (order 3)
- **Solution**: Use linear (order 1) or accept artifacts for higher quality

## Tips

1. **Always backup original data** before resampling
2. **Check spacing first**: `nibabel.load('file.nii.gz').header['pixdim'][1:4]`
3. **For label maps**: Always use `--order 0` (nearest neighbor)
4. **For upsampling**: Use `--order 3` (cubic) for better quality
5. **For downsampling**: Use `--order 1` (linear) - good enough
6. **Verify output**: Use `--compare` flag or check manually
7. **Consider disk space**: Upsampling can create VERY large files
8. **Memory usage**: Large volumes may need chunked processing (not in this script)

## Advanced: Understanding Zoom Factors

```
Zoom Factor = Original Spacing / Target Spacing

Examples:
100μm → 150μm: zoom = 100/150 = 0.667 (downsampling)
100μm → 50μm:  zoom = 100/50  = 2.0   (upsampling)
100μm → 100μm: zoom = 100/100 = 1.0   (no change)
```

The script calculates this automatically, but understanding helps predict output size.
