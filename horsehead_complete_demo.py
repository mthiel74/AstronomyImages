#!/usr/bin/env python3
"""
Complete Remote Telescope/Virtual Observatory API Guide
Demonstrates the full workflow with simulated data (network blocked in this environment)

WORKING APIS FOR REMOTE ASTRONOMY:
===================================

1. NASA SkyView Virtual Observatory
   URL: https://skyview.gsfc.nasa.gov/current/cgi/titlepage.pl
   API Endpoint: https://skyview.gsfc.nasa.gov/cgi-bin/images
   Authentication: None required
   Rate limit: Reasonable use
   
2. SDSS SkyServer
   URL: http://skyserver.sdss.org/
   API Endpoint: http://skyserver.sdss.org/dr18/SkyServerWS/ImgCutout/getjpeg
   Authentication: None for public data
   
3. Las Cumbres Observatory (LCO) - REAL TELESCOPE CONTROL
   URL: https://observe.lco.global/
   API: https://developers.lco.global/
   Authentication: API token required (free for education)
   Python library: pip install lcogt-comm
   
4. ESO Archive (European Southern Observatory)
   URL: http://archive.eso.org/
   API: Programmatic access via astroquery.eso
   
5. STScI (Space Telescope Science Institute)
   URL: https://archive.stsci.edu/
   API: astroquery.mast for Hubble/Webb data
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from astropy.io import fits
from astropy.wcs import WCS
import json

# Horsehead Nebula physical parameters
HORSEHEAD_RA = 85.2479   # degrees (5h 40m 59.5s)
HORSEHEAD_DEC = -2.4583  # degrees (-2° 27' 30")
HORSEHEAD_SIZE_ARCMIN = 8  # approximately 8 arcminutes tall

def create_simulated_horsehead_image(size=1000, noise_level=0.02):
    """
    Create simulated Horsehead Nebula observation
    Models the dark nebula as absorption feature against emission background
    
    Physics:
    - Background emission from IC 434 (H-alpha)
    - Foreground dust absorption (tau ~ 1-5 optical depths)
    - Scattered light and PSF effects
    """
    print("Simulating Horsehead Nebula observation...")
    print(f"  Image size: {size}x{size} pixels")
    print(f"  Noise level: {noise_level}")
    
    # Create coordinate grid
    x = np.linspace(-1, 1, size)
    y = np.linspace(-1, 1, size)
    X, Y = np.meshgrid(x, y)
    
    # Background emission nebula (IC 434) - H-alpha emission
    # Model as Gaussian cloud + structured filaments
    background = 100 + 150 * np.exp(-(X**2 + Y**2) / 0.5)
    
    # Add filamentary structure
    for i in range(5):
        angle = np.random.rand() * np.pi
        freq = np.random.rand() * 10 + 5
        amplitude = np.random.rand() * 30 + 10
        background += amplitude * np.sin(freq * (X * np.cos(angle) + Y * np.sin(angle)))
    
    # Create Horsehead shape (simplified)
    # The actual nebula resembles a horse's head in silhouette
    horsehead_mask = np.zeros((size, size))
    
    # "Head" region (absorption)
    head_x, head_y = 0.1, -0.2
    head_width, head_height = 0.3, 0.4
    horsehead_mask += np.exp(-((X - head_x)**2 / head_width**2 + (Y - head_y)**2 / head_height**2) * 8)
    
    # "Neck" region
    neck_x, neck_y = 0.15, 0.15
    neck_width, neck_height = 0.15, 0.35
    horsehead_mask += 0.7 * np.exp(-((X - neck_x)**2 / neck_width**2 + (Y - neck_y)**2 / neck_height**2) * 6)
    
    # Apply absorption: I_observed = I_background * exp(-tau)
    # tau is the optical depth through the dust
    optical_depth = horsehead_mask * 3.0  # Max optical depth ~3
    transmitted_fraction = np.exp(-optical_depth)
    
    # Final observed intensity
    image = background * transmitted_fraction
    
    # Add photon noise (Poisson statistics)
    image = np.random.poisson(image).astype(float)
    
    # Add readout noise (Gaussian)
    readout_noise = np.random.normal(0, noise_level * np.mean(image), image.shape)
    image += readout_noise
    
    # Ensure positive values
    image = np.maximum(image, 0)
    
    print(f"  Intensity range: [{image.min():.1f}, {image.max():.1f}] ADU")
    print(f"  Mean background: {background.mean():.1f} ADU")
    print(f"  Contrast ratio: {image.max()/image.min():.2f}")
    
    return image, background, optical_depth

def create_fits_with_wcs(image_data, ra, dec, pixel_scale=1.0):
    """
    Create FITS file with proper World Coordinate System (WCS)
    
    Parameters:
    -----------
    image_data : ndarray
        2D image array
    ra, dec : float
        Center coordinates in degrees
    pixel_scale : float
        Pixel scale in arcsec/pixel
    """
    # Create FITS HDU
    hdu = fits.PrimaryHDU(image_data)
    
    # Add WCS information to header
    header = hdu.header
    
    # Basic metadata
    header['OBJECT'] = 'Horsehead Nebula'
    header['TELESCOP'] = 'Simulated Observatory'
    header['INSTRUME'] = 'Virtual CCD'
    header['FILTER'] = 'H-alpha'
    header['EXPTIME'] = 300.0  # seconds
    header['DATE-OBS'] = '2025-11-17T00:00:00'
    
    # WCS keywords (FITS standard)
    ny, nx = image_data.shape
    header['CTYPE1'] = 'RA---TAN'  # Tangent plane projection
    header['CTYPE2'] = 'DEC--TAN'
    header['CRVAL1'] = ra   # Reference RA
    header['CRVAL2'] = dec  # Reference Dec
    header['CRPIX1'] = nx / 2.0  # Reference pixel X
    header['CRPIX2'] = ny / 2.0  # Reference pixel Y
    header['CDELT1'] = -pixel_scale / 3600.0  # degrees per pixel (negative for RA)
    header['CDELT2'] = pixel_scale / 3600.0   # degrees per pixel
    header['CUNIT1'] = 'deg'
    header['CUNIT2'] = 'deg'
    header['EQUINOX'] = 2000.0
    
    return hdu

def analyze_dark_nebula(image, threshold_percentile=30):
    """
    Quantitative analysis of dark nebula features
    
    Returns:
    --------
    stats : dict
        Statistical measurements
    """
    print("\n--- Dark Nebula Analysis ---")
    
    # Identify dark regions (absorption features)
    threshold = np.percentile(image, threshold_percentile)
    dark_mask = image < threshold
    
    # Calculate statistics
    dark_pixels = np.sum(dark_mask)
    total_pixels = image.size
    dark_fraction = dark_pixels / total_pixels
    
    background_mean = np.mean(image[~dark_mask])
    nebula_mean = np.mean(image[dark_mask])
    
    # Optical depth estimation: tau = -ln(I_observed / I_background)
    # Avoid log(0) by adding small epsilon
    intensity_ratio = (nebula_mean + 1e-6) / (background_mean + 1e-6)
    estimated_tau = -np.log(intensity_ratio)
    
    stats = {
        'dark_region_fraction': dark_fraction,
        'background_intensity': background_mean,
        'nebula_intensity': nebula_mean,
        'contrast': background_mean / nebula_mean,
        'estimated_optical_depth': estimated_tau,
        'extinction_magnitudes': estimated_tau * 1.086,  # Convert tau to magnitudes
    }
    
    print(f"Dark region coverage: {dark_fraction*100:.1f}% of image")
    print(f"Background intensity: {background_mean:.1f} ADU")
    print(f"Nebula intensity: {nebula_mean:.1f} ADU")
    print(f"Intensity contrast: {stats['contrast']:.2f}x")
    print(f"Estimated optical depth (tau): {estimated_tau:.2f}")
    print(f"Estimated extinction: {stats['extinction_magnitudes']:.2f} magnitudes")
    
    return stats, dark_mask

def plot_comprehensive_analysis(image, background, optical_depth, dark_mask):
    """Create comprehensive visualization"""
    
    fig = plt.figure(figsize=(16, 12))
    gs = fig.add_gridspec(3, 3, hspace=0.3, wspace=0.3)
    
    # 1. Original simulated observation
    ax1 = fig.add_subplot(gs[0, 0])
    im1 = ax1.imshow(image, cmap='gray', origin='lower')
    ax1.set_title('Simulated Observation\n(H-alpha filter)', fontsize=11)
    ax1.set_xlabel('Pixels (X)')
    ax1.set_ylabel('Pixels (Y)')
    plt.colorbar(im1, ax=ax1, label='Intensity (ADU)')
    
    # 2. Background emission (IC 434)
    ax2 = fig.add_subplot(gs[0, 1])
    im2 = ax2.imshow(background, cmap='hot', origin='lower')
    ax2.set_title('Background Emission Nebula\n(IC 434, intrinsic)', fontsize=11)
    ax2.set_xlabel('Pixels (X)')
    plt.colorbar(im2, ax=ax2, label='Intensity (ADU)')
    
    # 3. Optical depth map
    ax3 = fig.add_subplot(gs[0, 2])
    im3 = ax3.imshow(optical_depth, cmap='viridis', origin='lower')
    ax3.set_title('Dust Optical Depth\n(τ, absorption)', fontsize=11)
    ax3.set_xlabel('Pixels (X)')
    plt.colorbar(im3, ax=ax3, label='Optical Depth (τ)')
    
    # 4. Enhanced contrast view
    ax4 = fig.add_subplot(gs[1, 0])
    # Logarithmic scaling for detail
    log_image = np.log10(image + 1)
    im4 = ax4.imshow(log_image, cmap='gray', origin='lower')
    ax4.set_title('Log-Scaled Observation\n(enhanced contrast)', fontsize=11)
    ax4.set_xlabel('Pixels (X)')
    ax4.set_ylabel('Pixels (Y)')
    plt.colorbar(im4, ax=ax4, label='log₁₀(Intensity)')
    
    # 5. Segmented dark regions
    ax5 = fig.add_subplot(gs[1, 1])
    masked_view = np.ma.masked_where(~dark_mask, image)
    ax5.imshow(image, cmap='gray', origin='lower', alpha=0.5)
    im5 = ax5.imshow(masked_view, cmap='Reds', origin='lower', alpha=0.8)
    ax5.set_title('Segmented Dark Nebula\n(absorption regions)', fontsize=11)
    ax5.set_xlabel('Pixels (X)')
    plt.colorbar(im5, ax=ax5, label='Intensity (ADU)')
    
    # 6. Intensity histogram
    ax6 = fig.add_subplot(gs[1, 2])
    ax6.hist(image.flatten(), bins=100, color='steelblue', alpha=0.7, edgecolor='black')
    ax6.axvline(np.percentile(image, 30), color='red', linestyle='--', 
                label=f'30th percentile: {np.percentile(image, 30):.1f}')
    ax6.set_xlabel('Pixel Intensity (ADU)', fontsize=10)
    ax6.set_ylabel('Frequency', fontsize=10)
    ax6.set_title('Intensity Distribution', fontsize=11)
    ax6.legend()
    ax6.grid(True, alpha=0.3)
    
    # 7. Horizontal cross-section
    ax7 = fig.add_subplot(gs[2, 0])
    mid_row = image.shape[0] // 2
    ax7.plot(image[mid_row, :], 'b-', linewidth=1, label='Observed')
    ax7.plot(background[mid_row, :], 'r--', linewidth=1, alpha=0.7, label='Background')
    ax7.set_xlabel('Pixel Position (X)', fontsize=10)
    ax7.set_ylabel('Intensity (ADU)', fontsize=10)
    ax7.set_title(f'Horizontal Cross-Section (Y={mid_row})', fontsize=11)
    ax7.legend()
    ax7.grid(True, alpha=0.3)
    
    # 8. Vertical cross-section
    ax8 = fig.add_subplot(gs[2, 1])
    mid_col = image.shape[1] // 2
    ax8.plot(image[:, mid_col], 'b-', linewidth=1, label='Observed')
    ax8.plot(background[:, mid_col], 'r--', linewidth=1, alpha=0.7, label='Background')
    ax8.set_xlabel('Pixel Position (Y)', fontsize=10)
    ax8.set_ylabel('Intensity (ADU)', fontsize=10)
    ax8.set_title(f'Vertical Cross-Section (X={mid_col})', fontsize=11)
    ax8.legend()
    ax8.grid(True, alpha=0.3)
    
    # 9. 3D surface plot
    ax9 = fig.add_subplot(gs[2, 2], projection='3d')
    # Downsample for 3D visualization
    step = 20
    X = np.arange(0, image.shape[1], step)
    Y = np.arange(0, image.shape[0], step)
    X, Y = np.meshgrid(X, Y)
    Z = image[::step, ::step]
    surf = ax9.plot_surface(X, Y, Z, cmap='terrain', alpha=0.8, 
                            linewidth=0, antialiased=True)
    ax9.set_xlabel('X (pixels)', fontsize=8)
    ax9.set_ylabel('Y (pixels)', fontsize=8)
    ax9.set_zlabel('Intensity', fontsize=8)
    ax9.set_title('3D Intensity Surface', fontsize=11)
    ax9.view_init(elev=30, azim=45)
    
    fig.suptitle('Horsehead Nebula (Barnard 33) - Complete Analysis\n' + 
                 'Simulated H-alpha Observation', 
                 fontsize=14, fontweight='bold', y=0.995)
    
    return fig

def demonstrate_api_calls():
    """
    Show example API call structures for real services
    (These would work with proper network access)
    """
    
    print("\n" + "="*70)
    print("EXAMPLE API CALLS FOR REAL SERVICES")
    print("="*70)
    
    examples = {
        'SkyView': {
            'url': 'https://skyview.gsfc.nasa.gov/cgi-bin/images',
            'method': 'GET',
            'params': {
                'Position': f'{HORSEHEAD_RA},{HORSEHEAD_DEC}',
                'Survey': 'DSS2 Red',
                'Pixels': '1000,1000',
                'Size': 0.5,  # degrees
                'Return': 'FITS'
            },
            'python': """
import requests
response = requests.get(
    'https://skyview.gsfc.nasa.gov/cgi-bin/images',
    params={
        'Position': '85.2479,-2.4583',
        'Survey': 'DSS2 Red',
        'Pixels': '1000,1000',
        'Size': 0.5,
        'Return': 'FITS'
    }
)
fits_data = response.content
"""
        },
        'Las Cumbres Observatory': {
            'url': 'https://observe.lco.global/api/requestgroups/',
            'method': 'POST',
            'authentication': 'Token YOUR_API_TOKEN',
            'body': {
                'name': 'Horsehead Nebula Observation',
                'proposal': 'EDUCATION_PROJECT',
                'ipp_value': 1.0,
                'operator': 'SINGLE',
                'observation_type': 'NORMAL',
                'requests': [{
                    'target': {
                        'name': 'Horsehead Nebula',
                        'type': 'ICRS',
                        'ra': HORSEHEAD_RA,
                        'dec': HORSEHEAD_DEC
                    },
                    'configurations': [{
                        'type': 'EXPOSE',
                        'instrument_type': '1M0-SCICAM-SINISTRO',
                        'exposure_time': 300,
                        'exposure_count': 3,
                        'filter': 'rp'
                    }],
                    'windows': [{
                        'start': '2025-11-17T00:00:00Z',
                        'end': '2025-11-18T00:00:00Z'
                    }],
                    'location': {'telescope_class': '1m0'}
                }]
            },
            'python': """
import requests
headers = {'Authorization': 'Token YOUR_API_TOKEN'}
observation_request = {
    'name': 'Horsehead Nebula',
    'proposal': 'YOUR_PROPOSAL_ID',
    'ipp_value': 1.0,
    'operator': 'SINGLE',
    'observation_type': 'NORMAL',
    'requests': [{
        'target': {
            'name': 'Horsehead',
            'type': 'ICRS',
            'ra': 85.2479,
            'dec': -2.4583
        },
        'configurations': [{
            'type': 'EXPOSE',
            'instrument_type': '1M0-SCICAM-SINISTRO',
            'exposure_time': 300,
            'exposure_count': 3,
            'filter': 'rp'
        }],
        'windows': [{
            'start': '2025-11-17T00:00:00Z',
            'end': '2025-11-18T00:00:00Z'
        }],
        'location': {'telescope_class': '1m0'}
    }]
}
response = requests.post(
    'https://observe.lco.global/api/requestgroups/',
    headers=headers,
    json=observation_request
)
"""
        }
    }
    
    for service, details in examples.items():
        print(f"\n--- {service} ---")
        print(f"Endpoint: {details['url']}")
        print(f"Method: {details['method']}")
        if 'authentication' in details:
            print(f"Auth: {details['authentication']}")
        print(f"\nPython Example:")
        print(details['python'])
    
    return examples

def main():
    """Main demonstration"""
    
    print("="*70)
    print("REMOTE TELESCOPE & VIRTUAL OBSERVATORY DEMONSTRATION")
    print("="*70)
    print(f"\nTarget: Horsehead Nebula (Barnard 33)")
    print(f"Coordinates: RA={HORSEHEAD_RA}° ({85.2479/15:.2f}h), Dec={HORSEHEAD_DEC}°")
    print(f"Constellation: Orion")
    print(f"Type: Dark nebula (dust absorption feature)")
    print(f"Size: ~{HORSEHEAD_SIZE_ARCMIN}' x 5'")
    
    # Generate simulated observation
    print("\n" + "="*70)
    image, background, optical_depth = create_simulated_horsehead_image(size=800)
    
    # Save as FITS
    print("\n--- Creating FITS File ---")
    hdu = create_fits_with_wcs(image, HORSEHEAD_RA, HORSEHEAD_DEC, pixel_scale=1.0)
    fits_path = '/mnt/user-data/outputs/horsehead_simulated.fits'
    hdu.writeto(fits_path, overwrite=True)
    print(f"✓ FITS saved: {fits_path}")
    
    # Analyze
    stats, dark_mask = analyze_dark_nebula(image)
    
    # Save statistics
    stats_path = '/mnt/user-data/outputs/horsehead_analysis.json'
    with open(stats_path, 'w') as f:
        json.dump(stats, f, indent=2)
    print(f"✓ Statistics saved: {stats_path}")
    
    # Create visualizations
    print("\n--- Creating Visualizations ---")
    fig = plot_comprehensive_analysis(image, background, optical_depth, dark_mask)
    plot_path = '/mnt/user-data/outputs/horsehead_analysis.png'
    plt.savefig(plot_path, dpi=150, bbox_inches='tight')
    print(f"✓ Analysis plot saved: {plot_path}")
    plt.close()
    
    # Show API examples
    api_examples = demonstrate_api_calls()
    
    # Save API documentation
    api_doc_path = '/mnt/user-data/outputs/api_examples.json'
    with open(api_doc_path, 'w') as f:
        json.dump(api_examples, f, indent=2)
    print(f"\n✓ API documentation saved: {api_doc_path}")
    
    print("\n" + "="*70)
    print("DEMONSTRATION COMPLETE")
    print("="*70)
    print("\nFiles created:")
    print(f"  1. {fits_path} - FITS image with WCS")
    print(f"  2. {plot_path} - Comprehensive analysis")
    print(f"  3. {stats_path} - Quantitative measurements")
    print(f"  4. {api_doc_path} - API call examples")
    
    print("\nNext steps to use REAL telescopes:")
    print("  1. SkyView (free): Use code from horsehead_skyview.py")
    print("  2. LCO (free for education): Sign up at observe.lco.global")
    print("  3. iTelescope (paid): Professional-grade remote access")

if __name__ == "__main__":
    main()
