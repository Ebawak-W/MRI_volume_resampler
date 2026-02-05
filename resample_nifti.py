"""
NIfTI Volume Resampling Script
Resample a 100 micron or any isotropic volume to a new micron spacing
"""

import numpy as np
import nibabel as nib
from scipy.ndimage import zoom
import argparse


def resample_nifti_spacing(input_path, output_path, new_spacing=None,
                           new_x_spacing=None, new_y_spacing=None, new_z_spacing=None,
                           interpolation_order=1, verbose=True):
    """
    Resample a NIfTI volume by changing the voxel spacing in any direction(s).
    
    Parameters:
    -----------
    input_path : str
        Path to input NIfTI file
    output_path : str
        Path to save resampled NIfTI file
    new_spacing : tuple of float, optional
        New spacing as (x, y, z) in mm. If provided, overrides individual spacing params.
        Example: (0.1, 0.1, 0.15) for 100x100x150 microns
    new_x_spacing : float, optional
        New x-spacing in mm. If None, keeps original spacing.
    new_y_spacing : float, optional
        New y-spacing in mm. If None, keeps original spacing.
    new_z_spacing : float, optional
        New z-spacing in mm. If None, keeps original spacing.
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
        
    Examples:
    ---------
    # Using tuple
    resample_nifti_spacing('in.nii', 'out.nii', new_spacing=(0.1, 0.1, 0.15))
    
    # Using individual parameters (change only z)
    resample_nifti_spacing('in.nii', 'out.nii', new_z_spacing=0.15)
    
    # Change multiple directions
    resample_nifti_spacing('in.nii', 'out.nii', new_x_spacing=0.2, new_z_spacing=0.15)
    """
    
    # Load the NIfTI file
    if verbose:
        print(f"Loading NIfTI file: {input_path}")
    
    nifti_img = nib.load(input_path)
    data = nifti_img.get_fdata()
    affine = nifti_img.affine.copy()
    header = nifti_img.header.copy()
    
    # Get original spacing from the affine matrix
    # Spacing is the magnitude of each column vector (handles rotation)
    original_spacing = np.sqrt(np.sum(affine[:3, :3] ** 2, axis=0))
    
    if verbose:
        print(f"\nOriginal volume shape: {data.shape}")
        print(f"Original spacing (mm): x={original_spacing[0]:.4f}, "
              f"y={original_spacing[1]:.4f}, z={original_spacing[2]:.4f}")
        print(f"Original spacing (μm): x={original_spacing[0]*1000:.2f}, "
              f"y={original_spacing[1]*1000:.2f}, z={original_spacing[2]*1000:.2f}")
    
    # Determine target spacing for each axis
    if new_spacing is not None:
        # Use tuple input
        if len(new_spacing) != 3:
            raise ValueError("new_spacing must be a tuple of 3 values (x, y, z)")
        target_spacing = np.array(new_spacing, dtype=float)
        if verbose:
            print(f"\nUsing tuple input: new_spacing={new_spacing}")
    else:
        # Use individual parameters, defaulting to original spacing
        target_spacing = np.array([
            new_x_spacing if new_x_spacing is not None else original_spacing[0],
            new_y_spacing if new_y_spacing is not None else original_spacing[1],
            new_z_spacing if new_z_spacing is not None else original_spacing[2]  
        ], dtype=float)
        
        if verbose:
            changed_axes = []
            if new_x_spacing is not None:
                changed_axes.append(f"x={new_x_spacing}")
            if new_y_spacing is not None:
                changed_axes.append(f"y={new_y_spacing}")
            if new_z_spacing is not None:
                changed_axes.append(f"z={new_z_spacing}")
            # elif new_x_spacing is None and new_y_spacing is None:
            #     changed_axes.append(f"z=0.15 (default)")
            
            if changed_axes:
                print(f"\nChanging spacing: {', '.join(changed_axes)}")
    
    if verbose:
        print(f"\nTarget spacing (mm): x={target_spacing[0]:.4f}, "
              f"y={target_spacing[1]:.4f}, z={target_spacing[2]:.4f}")
        print(f"Target spacing (μm): x={target_spacing[0]*1000:.2f}, "
              f"y={target_spacing[1]*1000:.2f}, z={target_spacing[2]*1000:.2f}")
    
    # Calculate zoom factors (scaling factors)
    # zoom_factor = original_spacing / new_spacing
    zoom_factors = original_spacing / target_spacing
    
    if verbose:
        print(f"\nZoom factors: x={zoom_factors[0]:.4f}, y={zoom_factors[1]:.4f}, "
              f"z={zoom_factors[2]:.4f}")
        
        # Identify if upsampling or downsampling
        for i, (axis, zoom_val) in enumerate(zip(['x', 'y', 'z'], zoom_factors)):
            if zoom_val > 1.0:
                print(f"  {axis}-axis: Downsampling ({zoom_val:.4f}x)")
            elif zoom_val < 1.0:
                print(f"  {axis}-axis: Upsampling ({zoom_val:.4f}x)")
            else:
                print(f"  {axis}-axis: No change")
        
        interpolation_names = {
            0: "Nearest Neighbor",
            1: "Linear (Trilinear)",
            2: "Quadratic",
            3: "Cubic"
        }
        print(f"\nInterpolation method: {interpolation_names.get(interpolation_order, 'Unknown')}")
    
    # Perform resampling
    if verbose:
        print("\nResampling volume...")
    
    resampled_data = zoom(data, zoom_factors, order=interpolation_order, mode='nearest')
    
    if verbose:
        print(f"Resampled volume shape: {resampled_data.shape}")
        for i, axis in enumerate(['X', 'Y', 'Z']):
            old_size = data.shape[i]
            new_size = resampled_data.shape[i]
            change_pct = 100 * (new_size / old_size - 1)
            print(f"  {axis}: {old_size} → {new_size} ({change_pct:+.1f}%)")
    
    # Update the affine matrix for new spacing
    new_affine = affine.copy()
    
    # For rotated volumes, we need to scale each column vector
    for i in range(3):
        # Get the column vector
        col_vector = affine[:3, i]
        # Calculate its current magnitude (original spacing)
        current_magnitude = np.sqrt(np.sum(col_vector ** 2))
        # Scale it to have the new magnitude (target spacing)
        if current_magnitude > 0:
            scale_factor = target_spacing[i] / current_magnitude
            new_affine[:3, i] = col_vector * scale_factor
    
    # Update the header with new dimensions
    new_header = header.copy()
    new_header.set_data_shape(resampled_data.shape)
    
    # Update pixdim in header (voxel dimensions)
    pixdim = new_header['pixdim']
    pixdim[1:4] = target_spacing  # pixdim[1]=x, pixdim[2]=y, pixdim[3]=z
    new_header['pixdim'] = pixdim
    
    if verbose:
        # Verify new spacing
        verify_spacing = np.sqrt(np.sum(new_affine[:3, :3] ** 2, axis=0))
        print(f"\nVerified new spacing (mm): x={verify_spacing[0]:.4f}, "
              f"y={verify_spacing[1]:.4f}, z={verify_spacing[2]:.4f}")
        print(f"Verified new spacing (μm): x={verify_spacing[0]*1000:.2f}, "
              f"y={verify_spacing[1]*1000:.2f}, z={verify_spacing[2]*1000:.2f}")
    
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
    

def main():
    """Command-line interface for the resampling script."""
    
    parser = argparse.ArgumentParser(
        description='Resample NIfTI MRI volume to new voxel spacing',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
            Examples:
            # Change only z-spacing to 150 microns (linear interpolation)
            python resample_nifti.py input.nii.gz output.nii.gz --z-spacing 0.15
            
            # Change multiple axes
            python resample_nifti.py input.nii.gz output.nii.gz --x-spacing 0.2 --z-spacing 0.15
            
            # Use tuple input for all three axes
            python resample_nifti.py input.nii.gz output.nii.gz --spacing 0.1 0.1 0.15
            
            # With cubic interpolation
            python resample_nifti.py input.nii.gz output.nii.gz --z-spacing 0.15 --order 3
            
            # With comparison statistics
            python resample_nifti.py input.nii.gz output.nii.gz --z-spacing 0.15 --compare
        """
    )
    
    parser.add_argument('input', help='Input NIfTI file path')
    parser.add_argument('output', help='Output NIfTI file path')
    
    # Spacing options
    spacing_group = parser.add_mutually_exclusive_group()
    spacing_group.add_argument('--spacing', nargs=3, type=float, metavar=('X', 'Y', 'Z'),
                              help='New spacing as three values in mm (e.g., 0.1 0.1 0.15)')
    spacing_group.add_argument('--xyz-spacing', dest='spacing_tuple', nargs=3, type=float,
                              metavar=('X', 'Y', 'Z'),
                              help='Alternative syntax for --spacing')
    
    parser.add_argument('--x-spacing', type=float, default=None,
                       help='New x-spacing in mm (keeps original if not specified)')
    parser.add_argument('--y-spacing', type=float, default=None,
                       help='New y-spacing in mm (keeps original if not specified)')
    parser.add_argument('--z-spacing', type=float, default=None,
                       help='New z-spacing in mm (keeps original if not specified)')
    
    parser.add_argument('--order', type=int, default=1, choices=[0, 1, 2, 3],
                       help='Interpolation order: 0=nearest, 1=linear (default), '
                            '2=quadratic, 3=cubic')
    parser.add_argument('--compare', action='store_true',
                       help='Compare original and resampled volumes')
    parser.add_argument('--quiet', action='store_true',
                       help='Suppress verbose output')
    
    args = parser.parse_args()
    
    # Determine spacing parameters
    if args.spacing is not None:
        new_spacing = tuple(args.spacing)
        new_x = new_y = new_z = None
    elif args.spacing_tuple is not None:
        new_spacing = tuple(args.spacing_tuple)
        new_x = new_y = new_z = None
    else:
        new_spacing = None
        new_x = args.x_spacing
        new_y = args.y_spacing
        new_z = args.z_spacing
    
    # Perform resampling
    resample_nifti_spacing(
        args.input,
        args.output,
        new_spacing=new_spacing,
        new_x_spacing=new_x,
        new_y_spacing=new_y,
        new_z_spacing=new_z,
        interpolation_order=args.order,
        verbose=not args.quiet
    )
    
    # Compare volumes if requested
    if args.compare:
        compare_volumes(args.input, args.output)


if __name__ == '__main__':
    main()
