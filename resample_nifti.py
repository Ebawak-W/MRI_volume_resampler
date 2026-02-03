"""
NIfTI Volume Resampling Script
Resample a 100 micron isotropic volume to 100x100x150 micron spacing
"""

import numpy as np
import nibabel as nib
from scipy.ndimage import zoom
import argparse


def resample_nifti_z_spacing(input_path, output_path, new_z_spacing=0.15, 
                              interpolation_order=1, verbose=True):
    """
    Resample a NIfTI volume by changing only the z-direction spacing.
    
    Parameters:
    -----------
    input_path : str
        Path to input NIfTI file
    output_path : str
        Path to save resampled NIfTI file
    new_z_spacing : float
        New z-spacing in mm (0.15 mm = 150 microns)
    interpolation_order : int
        Interpolation order for scipy.ndimage.zoom
        0 = nearest neighbor
        1 = linear (default)
        2 = quadratic
        3 = cubic
    verbose : bool
        Print information about the resampling
    
    Returns:
    --------
    nibabel.Nifti1Image
        Resampled NIfTI image
    """
    
    # Load the NIfTI file
    if verbose:
        print(f"Loading NIfTI file: {input_path}")
    
    nifti_img = nib.load(input_path)
    data = nifti_img.get_fdata()
    affine = nifti_img.affine.copy()
    header = nifti_img.header.copy()
    
    # Get original spacing from the affine matrix
    original_spacing = np.abs(affine.diagonal()[:3])
    
    if verbose:
        print(f"\nOriginal volume shape: {data.shape}")
        print(f"Original spacing (mm): x={original_spacing[0]:.3f}, "
              f"y={original_spacing[1]:.3f}, z={original_spacing[2]:.3f}")
        print(f"Original spacing (μm): x={original_spacing[0]*1000:.1f}, "
              f"y={original_spacing[1]*1000:.1f}, z={original_spacing[2]*1000:.1f}")
    
    # Calculate zoom factors (scaling factors)
    # zoom_factor = original_spacing / new_spacing
    zoom_factors = [
        1.0,  # x-direction: no change
        1.0,  # y-direction: no change
        original_spacing[2] / new_z_spacing  # z-direction: downsample
    ]
    
    if verbose:
        print(f"\nTarget z-spacing: {new_z_spacing*1000:.1f} μm ({new_z_spacing:.3f} mm)")
        print(f"Zoom factors: x={zoom_factors[0]:.3f}, y={zoom_factors[1]:.3f}, "
              f"z={zoom_factors[2]:.3f}")
        
        interpolation_names = {
            0: "Nearest Neighbor",
            1: "Linear (Trilinear)",
            2: "Quadratic",
            3: "Cubic"
        }
        print(f"Interpolation method: {interpolation_names.get(interpolation_order, 'Unknown')}")
    
    # Perform resampling
    if verbose:
        print("\nResampling volume...")
    
    resampled_data = zoom(data, zoom_factors, order=interpolation_order, mode='nearest')
    
    if verbose:
        print(f"Resampled volume shape: {resampled_data.shape}")
        print(f"Size reduction in z: {data.shape[2]} → {resampled_data.shape[2]} slices "
              f"({100*(1 - resampled_data.shape[2]/data.shape[2]):.1f}% reduction)")
    
    # Update the affine matrix for new z-spacing
    new_affine = affine.copy()
    # Update z-spacing in the affine matrix (3rd diagonal element)
    new_affine[2, 2] = new_z_spacing if affine[2, 2] > 0 else -new_z_spacing
    
    # Update the header with new dimensions
    new_header = header.copy()
    new_header.set_data_shape(resampled_data.shape)
    
    # Update pixdim in header (voxel dimensions)
    pixdim = new_header['pixdim']
    pixdim[3] = new_z_spacing  # z-spacing (pixdim[1]=x, pixdim[2]=y, pixdim[3]=z)
    new_header['pixdim'] = pixdim
    
    if verbose:
        new_spacing = np.abs(new_affine.diagonal()[:3])
        print(f"\nNew spacing (mm): x={new_spacing[0]:.3f}, "
              f"y={new_spacing[1]:.3f}, z={new_spacing[2]:.3f}")
        print(f"New spacing (μm): x={new_spacing[0]*1000:.1f}, "
              f"y={new_spacing[1]*1000:.1f}, z={new_spacing[2]*1000:.1f}")
    
    # Create new NIfTI image
    resampled_nifti = nib.Nifti1Image(resampled_data, new_affine, new_header)
    
    # Save the resampled volume
    if verbose:
        print(f"\nSaving resampled volume to: {output_path}")
    
    nib.save(resampled_nifti, output_path)
    
    if verbose:
        print("Done!")
    
    return resampled_nifti


def compare_volumes(original_path, resampled_path):
    """
    Compare original and resampled volumes and print statistics.
    
    Parameters:
    -----------
    original_path : str
        Path to original NIfTI file
    resampled_path : str
        Path to resampled NIfTI file
    """
    
    print("\n" + "="*60)
    print("VOLUME COMPARISON")
    print("="*60)
    
    # Load both volumes
    orig_img = nib.load(original_path)
    resamp_img = nib.load(resampled_path)
    
    orig_data = orig_img.get_fdata()
    resamp_data = resamp_img.get_fdata()
    
    # Basic statistics
    print(f"\nOriginal volume:")
    print(f"  Shape: {orig_data.shape}")
    print(f"  Total voxels: {orig_data.size:,}")
    print(f"  Intensity range: [{orig_data.min():.2f}, {orig_data.max():.2f}]")
    print(f"  Mean intensity: {orig_data.mean():.2f}")
    print(f"  Std intensity: {orig_data.std():.2f}")
    
    print(f"\nResampled volume:")
    print(f"  Shape: {resamp_data.shape}")
    print(f"  Total voxels: {resamp_data.size:,}")
    print(f"  Intensity range: [{resamp_data.min():.2f}, {resamp_data.max():.2f}]")
    print(f"  Mean intensity: {resamp_data.mean():.2f}")
    print(f"  Std intensity: {resamp_data.std():.2f}")
    
    print(f"\nVoxel count change: {100*(resamp_data.size/orig_data.size - 1):.1f}%")
    print(f"Memory footprint change: "
          f"{100*(resamp_data.nbytes/orig_data.nbytes - 1):.1f}%")


def main():
    """Command-line interface for the resampling script."""
    
    parser = argparse.ArgumentParser(
        description='Resample NIfTI MRI volume from isotropic to anisotropic spacing',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage (linear interpolation)
  python resample_nifti.py input.nii.gz output.nii.gz
  
  # With custom z-spacing (200 microns)
  python resample_nifti.py input.nii.gz output.nii.gz --z-spacing 0.20
  
  # Using cubic interpolation
  python resample_nifti.py input.nii.gz output.nii.gz --order 3
  
  # With comparison statistics
  python resample_nifti.py input.nii.gz output.nii.gz --compare
        """
    )
    
    parser.add_argument('input', help='Input NIfTI file path')
    parser.add_argument('output', help='Output NIfTI file path')
    parser.add_argument('--z-spacing', type=float, default=0.15,
                       help='New z-spacing in mm (default: 0.15 = 150 microns)')
    parser.add_argument('--order', type=int, default=1, choices=[0, 1, 2, 3],
                       help='Interpolation order: 0=nearest, 1=linear (default), '
                            '2=quadratic, 3=cubic')
    parser.add_argument('--compare', action='store_true',
                       help='Compare original and resampled volumes')
    parser.add_argument('--quiet', action='store_true',
                       help='Suppress verbose output')
    
    args = parser.parse_args()
    
    # Perform resampling
    resample_nifti_z_spacing(
        args.input,
        args.output,
        new_z_spacing=args.z_spacing,
        interpolation_order=args.order,
        verbose=not args.quiet
    )
    
    # Compare volumes if requested
    if args.compare:
        compare_volumes(args.input, args.output)


if __name__ == '__main__':
    main()
