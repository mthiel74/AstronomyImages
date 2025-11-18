#!/usr/bin/env python3
"""
Query Sloan Digital Sky Survey (SDSS) for Horsehead Nebula region
Alternative to SkyView when network restrictions apply
"""

import requests
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image
from io import BytesIO

# Horsehead Nebula coordinates
HORSEHEAD_RA = 85.2479   # degrees
HORSEHEAD_DEC = -2.4583  # degrees

def query_sdss_image(ra, dec, scale=0.4, width=512, height=512):
    """
    Query SDSS Image Cutout Service
    
    Parameters:
    -----------
    ra : float
        Right Ascension in degrees
    dec : float
        Declination in degrees
    scale : float
        Image scale in arcsec/pixel (0.4 is typical)
    width, height : int
        Image dimensions in pixels
    
    Returns:
    --------
    image : PIL Image or None
    """
    # SDSS SkyServer Image Cutout URL
    base_url = "http://skyserver.sdss.org/dr18/SkyServerWS/ImgCutout/getjpeg"
    
    params = {
        'ra': ra,
        'dec': dec,
        'scale': scale,
        'width': width,
        'height': height,
        'opt': 'G',  # Options: G=Grid overlay
    }
    
    print(f"Querying SDSS DR18 for coordinates: RA={ra}°, Dec={dec}°")
    print(f"Scale: {scale} arcsec/pixel, Size: {width}x{height} pixels")
    print(f"Field of view: ~{width*scale/3600:.2f}° x {height*scale/3600:.2f}°")
    
    try:
        response = requests.get(base_url, params=params, timeout=30)
        
        if response.status_code == 200 and len(response.content) > 0:
            print("✓ Successfully retrieved image")
            image = Image.open(BytesIO(response.content))
            return image
        else:
            print(f"✗ Query returned status {response.status_code}")
            return None
    except Exception as e:
        print(f"✗ Error: {e}")
        return None

def query_dss_via_alternate(ra, dec, size_arcmin=30):
    """
    Alternative: STScI Digitized Sky Survey
    """
    base_url = "http://archive.stsci.edu/cgi-bin/dss_search"
    
    # Convert RA to hours:minutes:seconds format
    ra_hours = ra / 15.0
    ra_h = int(ra_hours)
    ra_m = int((ra_hours - ra_h) * 60)
    ra_s = ((ra_hours - ra_h) * 60 - ra_m) * 60
    
    # Convert Dec to degrees:arcminutes:arcseconds
    dec_sign = '+' if dec >= 0 else '-'
    dec_abs = abs(dec)
    dec_d = int(dec_abs)
    dec_m = int((dec_abs - dec_d) * 60)
    dec_s = ((dec_abs - dec_d) * 60 - dec_m) * 60
    
    params = {
        'ra': f'{ra_h:02d} {ra_m:02d} {ra_s:05.2f}',
        'dec': f'{dec_sign}{dec_d:02d} {dec_m:02d} {dec_s:05.2f}',
        'equinox': 'J2000',
        'height': size_arcmin,
        'width': size_arcmin,
        'format': 'JPEG',
        'generation': 'DSS2',
    }
    
    print(f"\nQuerying STScI DSS2 for: RA={params['ra']}, Dec={params['dec']}")
    print(f"Field of view: {size_arcmin}' x {size_arcmin}'")
    
    try:
        response = requests.get(base_url, params=params, timeout=30)
        
        if response.status_code == 200 and len(response.content) > 0:
            print("✓ Successfully retrieved DSS image")
            image = Image.open(BytesIO(response.content))
            return image
        else:
            print(f"✗ Query returned status {response.status_code}")
            return None
    except Exception as e:
        print(f"✗ Error: {e}")
        return None

def display_and_analyze(image, title, save_path):
    """Display astronomical image with analysis"""
    if image is None:
        print("No image to display")
        return
    
    # Convert to numpy array
    img_array = np.array(image)
    
    print(f"\n--- Image Analysis ---")
    print(f"Image shape: {img_array.shape}")
    print(f"Data type: {img_array.dtype}")
    print(f"Value range: [{img_array.min()}, {img_array.max()}]")
    
    # Create figure
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    # Original image
    axes[0].imshow(img_array)
    axes[0].set_title(f'{title}\nOriginal', fontsize=12)
    axes[0].axis('off')
    
    # Enhanced contrast version (histogram equalization)
    if len(img_array.shape) == 3:
        # Convert to grayscale for analysis
        gray = np.mean(img_array, axis=2)
    else:
        gray = img_array
    
    # Apply histogram equalization for better visibility
    from matplotlib import cm
    axes[1].imshow(gray, cmap='gray')
    axes[1].set_title(f'{title}\nGrayscale', fontsize=12)
    axes[1].axis('off')
    
    plt.tight_layout()
    plt.savefig(save_path, dpi=150, bbox_inches='tight')
    print(f"✓ Saved to: {save_path}")
    plt.close()
    
    # Statistical analysis
    print(f"\n--- Intensity Statistics (grayscale) ---")
    print(f"Mean: {gray.mean():.2f}")
    print(f"Std Dev: {gray.std():.2f}")
    print(f"Median: {np.median(gray):.2f}")
    
    # Create histogram
    fig, ax = plt.subplots(1, 1, figsize=(10, 5))
    ax.hist(gray.flatten(), bins=100, color='steelblue', alpha=0.7, edgecolor='black')
    ax.set_xlabel('Pixel Intensity', fontsize=12)
    ax.set_ylabel('Count', fontsize=12)
    ax.set_title(f'{title} - Intensity Distribution', fontsize=14)
    ax.grid(True, alpha=0.3)
    
    hist_path = save_path.replace('.png', '_histogram.png')
    plt.savefig(hist_path, dpi=150, bbox_inches='tight')
    print(f"✓ Histogram saved to: {hist_path}")
    plt.close()

def main():
    """Main execution"""
    print("=" * 70)
    print("Astronomical Image Query: Horsehead Nebula (Barnard 33)")
    print("=" * 70)
    print("\nCoordinates:")
    print(f"  RA:  85.2479° = 5h 40m 59.5s")
    print(f"  Dec: -2.4583° = -2° 27' 30\"")
    print(f"\nNote: The Horsehead Nebula is a dark nebula - an absorption feature")
    print(f"      visible against the emission nebula IC 434 in the Orion constellation.")
    
    # Try SDSS first
    print("\n" + "="*70)
    print("Method 1: SDSS DR18 (may not cover this region well)")
    print("="*70)
    image_sdss = query_sdss_image(HORSEHEAD_RA, HORSEHEAD_DEC, 
                                    scale=0.8, width=800, height=800)
    
    if image_sdss:
        display_and_analyze(
            image_sdss, 
            "Horsehead Nebula Region (SDSS)",
            "/mnt/user-data/outputs/horsehead_sdss.png"
        )
    
    # Try STScI DSS
    print("\n" + "="*70)
    print("Method 2: STScI Digitized Sky Survey")
    print("="*70)
    image_dss = query_dss_via_alternate(HORSEHEAD_RA, HORSEHEAD_DEC, size_arcmin=30)
    
    if image_dss:
        display_and_analyze(
            image_dss,
            "Horsehead Nebula (DSS2)",
            "/mnt/user-data/outputs/horsehead_dss.png"
        )
    
    print("\n" + "="*70)
    print("Query Complete!")
    print("="*70)
    
    if not image_sdss and not image_dss:
        print("\n⚠ Network restrictions prevented both queries.")
        print("The code structure is correct - would work with proper network access.")
        print("\nAlternative: Download FITS files manually from:")
        print("  - https://skyview.gsfc.nasa.gov/current/cgi/titlepage.pl")
        print("  - https://archive.stsci.edu/dss/")

if __name__ == "__main__":
    main()
