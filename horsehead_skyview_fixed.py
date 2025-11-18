#!/usr/bin/env python3
"""
Query NASA SkyView Virtual Observatory for Horsehead Nebula
Fixed version for local filesystem
"""

import requests
from astropy.io import fits
from astropy.wcs import WCS
import matplotlib.pyplot as plt
import numpy as np
from io import BytesIO
import os

# Horsehead Nebula coordinates
HORSEHEAD_RA = 85.2479
HORSEHEAD_DEC = -2.4583

def query_skyview(ra, dec, survey='DSS', pixels=800, fov=0.5):
    """Query SkyView for astronomical image"""
    base_url = "https://skyview.gsfc.nasa.gov/cgi-bin/images"
    
    params = {
        'Position': f'{ra},{dec}',
        'Survey': survey,
        'Pixels': f'{pixels},{pixels}',
        'Size': fov,
        'Return': 'FITS',
        'Scaling': 'Linear',
    }
    
    print(f"Querying SkyView for coordinates: RA={ra}°, Dec={dec}°")
    print(f"Survey: {survey}, FOV: {fov}°, Resolution: {pixels}x{pixels} pixels")
    
    response = requests.get(base_url, params=params, timeout=60)
    
    if response.status_code == 200:
        print("✓ Successfully retrieved data")
        return response.content
    else:
        raise Exception(f"Query failed with status code {response.status_code}")

def display_fits_image(fits_data, title="Astronomical Image"):
    """Display FITS image with proper scaling"""
    with fits.open(BytesIO(fits_data)) as hdul:
        image_data = hdul[0].data
        header = hdul[0].header
        wcs = WCS(header)
        
        print("\n--- FITS Header Information ---")
        print(f"Image shape: {image_data.shape}")
        print(f"Min value: {np.nanmin(image_data):.2e}")
        print(f"Max value: {np.nanmax(image_data):.2e}")
        print(f"Mean value: {np.nanmean(image_data):.2e}")
        
        # Create figure
        fig = plt.figure(figsize=(12, 10))
        ax = fig.add_subplot(111, projection=wcs)
        
        # Apply percentile clipping
        vmin = np.nanpercentile(image_data, 1)
        vmax = np.nanpercentile(image_data, 99.5)
        
        im = ax.imshow(image_data, cmap='gray', origin='lower',
                      vmin=vmin, vmax=vmax, interpolation='nearest')
        
        ax.coords.grid(True, color='cyan', ls=':', alpha=0.5)
        ax.coords[0].set_axislabel('Right Ascension (J2000)')
        ax.coords[1].set_axislabel('Declination (J2000)')
        
        cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        cbar.set_label('Intensity (ADU)', rotation=270, labelpad=20)
        
        ax.set_title(title, fontsize=14, pad=20)
        plt.tight_layout()
        
        return fig, ax, image_data, header

def save_fits(fits_data, filename):
    """Save FITS data to file in current directory"""
    with open(filename, 'wb') as f:
        f.write(fits_data)
    print(f"✓ FITS file saved to: {filename}")

def main():
    """Main execution"""
    print("=" * 60)
    print("NASA SkyView Query: Horsehead Nebula (Barnard 33)")
    print("=" * 60)
    
    # Query different surveys
    surveys = [
        ('DSS2 Red', 'Digitized Sky Survey (Red)', 0.5),
        ('2MASS-J', '2MASS Near-Infrared (J-band)', 0.5),
    ]
    
    for survey_name, description, fov in surveys:
        print(f"\n{'='*60}")
        print(f"Querying: {description}")
        print(f"{'='*60}")
        
        try:
            # Query SkyView
            fits_data = query_skyview(
                HORSEHEAD_RA, 
                HORSEHEAD_DEC, 
                survey=survey_name,
                pixels=1000,
                fov=fov
            )
            
            # Save FITS file (current directory)
            safe_name = survey_name.replace(' ', '_').replace('-', '_')
            filename = f'horsehead_{safe_name}.fits'
            save_fits(fits_data, filename)
            
            # Display image
            title = f"Horsehead Nebula - {description}"
            fig, ax, data, header = display_fits_image(fits_data, title)
            
            # Save plot (current directory)
            plot_filename = f'horsehead_{safe_name}.png'
            plt.savefig(plot_filename, dpi=150, bbox_inches='tight')
            print(f"✓ Plot saved to: {plot_filename}")
            plt.close()
            
        except Exception as e:
            print(f"✗ Error querying {survey_name}: {e}")
    
    print("\n" + "="*60)
    print("Analysis complete!")
    print("="*60)

if __name__ == "__main__":
    main()
