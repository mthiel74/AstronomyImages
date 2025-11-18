#!/usr/bin/env python3
"""
Interactive Astronomy Image Viewer
Click on sky map to download and view astronomical images
"""

import tkinter as tk
from tkinter import ttk, messagebox, filedialog
import threading
import requests
from astropy.io import fits
from astropy.wcs import WCS
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
import numpy as np
from io import BytesIO
import json

# Popular astronomical objects database
CATALOG = {
    'Horsehead Nebula': {'ra': 85.2479, 'dec': -2.4583, 'size': 0.3},
    'Orion Nebula (M42)': {'ra': 83.8221, 'dec': -5.3911, 'size': 1.0},
    'Andromeda Galaxy (M31)': {'ra': 10.6847, 'dec': 41.2687, 'size': 3.0},
    'Crab Nebula (M1)': {'ra': 83.6333, 'dec': 22.0145, 'size': 0.3},
    'Ring Nebula (M57)': {'ra': 283.3963, 'dec': 33.0295, 'size': 0.1},
    'Whirlpool Galaxy (M51)': {'ra': 202.4696, 'dec': 47.1952, 'size': 0.3},
    'Eagle Nebula (M16)': {'ra': 274.7000, 'dec': -13.8167, 'size': 0.5},
    'Lagoon Nebula (M8)': {'ra': 270.9167, 'dec': -24.3833, 'size': 1.0},
    'Pleiades (M45)': {'ra': 56.75, 'dec': 24.1167, 'size': 2.0},
    'Dumbbell Nebula (M27)': {'ra': 299.9013, 'dec': 22.7211, 'size': 0.3},
}

class AstronomyViewer:
    def __init__(self, root):
        self.root = root
        self.root.title("Interactive Astronomy Image Viewer")
        self.root.geometry("1400x900")
        
        self.current_fits_data = None
        self.current_image_data = None
        
        self.setup_ui()
        
    def setup_ui(self):
        """Create the user interface"""
        
        # Main container
        main_frame = ttk.Frame(self.root, padding="10")
        main_frame.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
        
        self.root.columnconfigure(0, weight=1)
        self.root.rowconfigure(0, weight=1)
        main_frame.columnconfigure(1, weight=1)
        main_frame.rowconfigure(1, weight=1)
        
        # Left panel - Controls
        control_frame = ttk.LabelFrame(main_frame, text="Controls", padding="10")
        control_frame.grid(row=0, column=0, rowspan=2, sticky=(tk.W, tk.E, tk.N, tk.S), padx=(0, 10))
        
        # Object catalog
        ttk.Label(control_frame, text="Popular Objects:").grid(row=0, column=0, sticky=tk.W, pady=5)
        self.object_var = tk.StringVar()
        object_dropdown = ttk.Combobox(control_frame, textvariable=self.object_var, 
                                       values=list(CATALOG.keys()), width=25, state='readonly')
        object_dropdown.grid(row=1, column=0, sticky=(tk.W, tk.E), pady=5)
        object_dropdown.current(0)
        
        ttk.Button(control_frame, text="Load Object", 
                  command=self.load_catalog_object).grid(row=2, column=0, sticky=(tk.W, tk.E), pady=5)
        
        # Separator
        ttk.Separator(control_frame, orient='horizontal').grid(row=3, column=0, sticky=(tk.W, tk.E), pady=10)
        
        # Manual coordinates
        ttk.Label(control_frame, text="Manual Coordinates:", font=('TkDefaultFont', 10, 'bold')).grid(row=4, column=0, sticky=tk.W, pady=5)
        
        ttk.Label(control_frame, text="Right Ascension (°):").grid(row=5, column=0, sticky=tk.W)
        self.ra_var = tk.StringVar(value="85.2479")
        ttk.Entry(control_frame, textvariable=self.ra_var, width=25).grid(row=6, column=0, sticky=(tk.W, tk.E), pady=2)
        
        ttk.Label(control_frame, text="Declination (°):").grid(row=7, column=0, sticky=tk.W)
        self.dec_var = tk.StringVar(value="-2.4583")
        ttk.Entry(control_frame, textvariable=self.dec_var, width=25).grid(row=8, column=0, sticky=(tk.W, tk.E), pady=2)
        
        ttk.Label(control_frame, text="Field of View (°):").grid(row=9, column=0, sticky=tk.W)
        self.fov_var = tk.StringVar(value="0.5")
        ttk.Entry(control_frame, textvariable=self.fov_var, width=25).grid(row=10, column=0, sticky=(tk.W, tk.E), pady=2)
        
        # Survey selection
        ttk.Label(control_frame, text="Survey:").grid(row=11, column=0, sticky=tk.W, pady=(10, 0))
        self.survey_var = tk.StringVar(value="DSS2 Red")
        surveys = ['DSS2 Red', 'DSS2 Blue', '2MASS-J', '2MASS-H', 'WISE 3.4']
        ttk.Combobox(control_frame, textvariable=self.survey_var, 
                    values=surveys, width=25, state='readonly').grid(row=12, column=0, sticky=(tk.W, tk.E), pady=2)
        
        # Resolution
        ttk.Label(control_frame, text="Resolution (pixels):").grid(row=13, column=0, sticky=tk.W)
        self.resolution_var = tk.StringVar(value="800")
        ttk.Combobox(control_frame, textvariable=self.resolution_var, 
                    values=['500', '800', '1000', '1500', '2000'], width=25, state='readonly').grid(row=14, column=0, sticky=(tk.W, tk.E), pady=2)
        
        # Download button
        ttk.Button(control_frame, text="📡 Download Image", 
                  command=self.download_image, 
                  style='Accent.TButton').grid(row=15, column=0, sticky=(tk.W, tk.E), pady=15)
        
        # Separator
        ttk.Separator(control_frame, orient='horizontal').grid(row=16, column=0, sticky=(tk.W, tk.E), pady=10)
        
        # Image controls
        ttk.Label(control_frame, text="Image Display:", font=('TkDefaultFont', 10, 'bold')).grid(row=17, column=0, sticky=tk.W, pady=5)
        
        ttk.Label(control_frame, text="Color Map:").grid(row=18, column=0, sticky=tk.W)
        self.cmap_var = tk.StringVar(value="gray")
        cmaps = ['gray', 'viridis', 'hot', 'cool', 'rainbow', 'plasma']
        ttk.Combobox(control_frame, textvariable=self.cmap_var, 
                    values=cmaps, width=25, state='readonly').grid(row=19, column=0, sticky=(tk.W, tk.E), pady=2)
        
        ttk.Button(control_frame, text="Update Display", 
                  command=self.update_display).grid(row=20, column=0, sticky=(tk.W, tk.E), pady=5)
        
        ttk.Button(control_frame, text="💾 Save FITS", 
                  command=self.save_fits).grid(row=21, column=0, sticky=(tk.W, tk.E), pady=2)
        
        ttk.Button(control_frame, text="💾 Save PNG", 
                  command=self.save_png).grid(row=22, column=0, sticky=(tk.W, tk.E), pady=2)
        
        # Status
        self.status_label = ttk.Label(control_frame, text="Ready", foreground="green")
        self.status_label.grid(row=23, column=0, sticky=tk.W, pady=10)
        
        # Right panel - Image display
        display_frame = ttk.LabelFrame(main_frame, text="Image View", padding="10")
        display_frame.grid(row=0, column=1, rowspan=2, sticky=(tk.W, tk.E, tk.N, tk.S))
        display_frame.columnconfigure(0, weight=1)
        display_frame.rowconfigure(0, weight=1)
        
        # Create matplotlib figure
        self.fig = Figure(figsize=(10, 8), dpi=100)
        self.ax = self.fig.add_subplot(111)
        self.canvas = FigureCanvasTkAgg(self.fig, master=display_frame)
        self.canvas.get_tk_widget().grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
        
        # Info panel at top
        info_frame = ttk.LabelFrame(main_frame, text="Image Information", padding="10")
        info_frame.grid(row=0, column=1, sticky=(tk.W, tk.E, tk.N), pady=(0, 10))
        
        self.info_text = tk.Text(info_frame, height=4, width=80, wrap=tk.WORD)
        self.info_text.grid(row=0, column=0, sticky=(tk.W, tk.E))
        self.info_text.insert('1.0', "No image loaded. Select an object or enter coordinates and click 'Download Image'.")
        self.info_text.config(state='disabled')
        
        # Initial display
        self.show_welcome_screen()
        
    def show_welcome_screen(self):
        """Display welcome message"""
        self.ax.clear()
        self.ax.text(0.5, 0.5, 'Interactive Astronomy Viewer\n\n' + 
                    'Select an object or enter coordinates\nthen click "Download Image"',
                    ha='center', va='center', fontsize=14, transform=self.ax.transAxes)
        self.ax.axis('off')
        self.canvas.draw()
        
    def load_catalog_object(self):
        """Load coordinates from catalog"""
        obj_name = self.object_var.get()
        if obj_name in CATALOG:
            obj = CATALOG[obj_name]
            self.ra_var.set(str(obj['ra']))
            self.dec_var.set(str(obj['dec']))
            self.fov_var.set(str(obj['size']))
            self.update_info(f"Loaded: {obj_name}\nRA: {obj['ra']}°, Dec: {obj['dec']}°, FOV: {obj['size']}°")
    
    def download_image(self):
        """Download image from SkyView"""
        try:
            ra = float(self.ra_var.get())
            dec = float(self.dec_var.get())
            fov = float(self.fov_var.get())
            survey = self.survey_var.get()
            resolution = int(self.resolution_var.get())
            
            self.status_label.config(text="Downloading...", foreground="orange")
            self.root.update()
            
            # Run in thread to avoid freezing UI
            thread = threading.Thread(target=self._download_thread, 
                                     args=(ra, dec, fov, survey, resolution))
            thread.daemon = True
            thread.start()
            
        except ValueError as e:
            messagebox.showerror("Invalid Input", f"Please check your coordinate values:\n{e}")
    
    def _download_thread(self, ra, dec, fov, survey, resolution):
        """Background download thread"""
        try:
            url = "https://skyview.gsfc.nasa.gov/cgi-bin/images"
            params = {
                'Position': f'{ra},{dec}',
                'Survey': survey,
                'Pixels': f'{resolution},{resolution}',
                'Size': fov,
                'Return': 'FITS',
                'Scaling': 'Linear',
            }
            
            response = requests.get(url, params=params, timeout=60)
            
            if response.status_code == 200 and len(response.content) > 1000:
                self.current_fits_data = response.content
                
                # Parse FITS
                with fits.open(BytesIO(self.current_fits_data)) as hdul:
                    self.current_image_data = hdul[0].data
                    header = hdul[0].header
                
                # Update UI in main thread
                self.root.after(0, self._display_image, header, survey, ra, dec)
                self.root.after(0, lambda: self.status_label.config(text="✓ Download complete", foreground="green"))
            else:
                self.root.after(0, lambda: messagebox.showerror("Download Failed", 
                               f"Failed to retrieve image. Status: {response.status_code}"))
                self.root.after(0, lambda: self.status_label.config(text="✗ Download failed", foreground="red"))
                
        except Exception as e:
            self.root.after(0, lambda: messagebox.showerror("Error", str(e)))
            self.root.after(0, lambda: self.status_label.config(text="✗ Error", foreground="red"))
    
    def _display_image(self, header, survey, ra, dec):
        """Display the downloaded image"""
        # Clear the entire figure to remove old colorbars
        self.fig.clear()
        self.ax = self.fig.add_subplot(111)
        
        # Apply colormap and scaling
        vmin = np.nanpercentile(self.current_image_data, 1)
        vmax = np.nanpercentile(self.current_image_data, 99)
        
        im = self.ax.imshow(self.current_image_data, 
                           cmap=self.cmap_var.get(), 
                           origin='lower',
                           vmin=vmin, 
                           vmax=vmax,
                           interpolation='nearest')
        
        self.fig.colorbar(im, ax=self.ax, label='Intensity (ADU)')
        
        # Add title
        self.ax.set_title(f'{survey} - RA: {ra:.4f}°, Dec: {dec:.4f}°', fontsize=12)
        self.ax.set_xlabel('X (pixels)')
        self.ax.set_ylabel('Y (pixels)')
        
        self.canvas.draw()
        
        # Update info
        info = f"Survey: {survey}\n"
        info += f"Coordinates: RA={ra:.4f}°, Dec={dec:.4f}°\n"
        info += f"Image size: {self.current_image_data.shape}\n"
        info += f"Intensity range: [{np.nanmin(self.current_image_data):.2e}, {np.nanmax(self.current_image_data):.2e}] ADU"
        self.update_info(info)
    
    def update_display(self):
        """Update display with new colormap"""
        if self.current_image_data is not None:
            ra = float(self.ra_var.get())
            dec = float(self.dec_var.get())
            survey = self.survey_var.get()
            
            # Clear the entire figure to remove old colorbars
            self.fig.clear()
            self.ax = self.fig.add_subplot(111)
            
            vmin = np.nanpercentile(self.current_image_data, 1)
            vmax = np.nanpercentile(self.current_image_data, 99)
            
            im = self.ax.imshow(self.current_image_data, 
                               cmap=self.cmap_var.get(), 
                               origin='lower',
                               vmin=vmin, 
                               vmax=vmax,
                               interpolation='nearest')
            
            self.fig.colorbar(im, ax=self.ax, label='Intensity (ADU)')
            self.ax.set_title(f'{survey} - RA: {ra:.4f}°, Dec: {dec:.4f}°', fontsize=12)
            self.ax.set_xlabel('X (pixels)')
            self.ax.set_ylabel('Y (pixels)')
            
            self.canvas.draw()
        else:
            messagebox.showinfo("No Image", "Please download an image first.")
    
    def save_fits(self):
        """Save FITS file"""
        if self.current_fits_data is None:
            messagebox.showinfo("No Image", "Please download an image first.")
            return
        
        filename = filedialog.asksaveasfilename(
            defaultextension=".fits",
            filetypes=[("FITS files", "*.fits"), ("All files", "*.*")],
            initialfile=f"astronomy_{self.survey_var.get().replace(' ', '_')}.fits"
        )
        
        if filename:
            with open(filename, 'wb') as f:
                f.write(self.current_fits_data)
            messagebox.showinfo("Saved", f"FITS file saved to:\n{filename}")
    
    def save_png(self):
        """Save PNG file"""
        if self.current_image_data is None:
            messagebox.showinfo("No Image", "Please download an image first.")
            return
        
        filename = filedialog.asksaveasfilename(
            defaultextension=".png",
            filetypes=[("PNG files", "*.png"), ("All files", "*.*")],
            initialfile=f"astronomy_{self.survey_var.get().replace(' ', '_')}.png"
        )
        
        if filename:
            self.fig.savefig(filename, dpi=150, bbox_inches='tight')
            messagebox.showinfo("Saved", f"PNG file saved to:\n{filename}")
    
    def update_info(self, text):
        """Update info panel"""
        self.info_text.config(state='normal')
        self.info_text.delete('1.0', tk.END)
        self.info_text.insert('1.0', text)
        self.info_text.config(state='disabled')

def main():
    root = tk.Tk()
    app = AstronomyViewer(root)
    root.mainloop()

if __name__ == "__main__":
    main()
