#!/usr/bin/env python3
"""
Interactive Sky Map Astronomy Viewer
Click anywhere on the sky to download and view that region
"""

import tkinter as tk
from tkinter import ttk, messagebox, filedialog
import threading
import requests
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import astropy.units as u
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure
import numpy as np
from io import BytesIO

class InteractiveSkyViewer:
    def __init__(self, root):
        self.root = root
        self.root.title("Interactive Sky Map - Click to Explore")
        self.root.geometry("1600x1000")
        
        self.current_fits_data = None
        self.current_image_data = None
        self.selected_ra = None
        self.selected_dec = None
        
        # Constellation data for sky map (simplified major constellations)
        self.setup_constellation_data()
        
        self.setup_ui()
        self.draw_sky_map()
        
    def setup_constellation_data(self):
        """Define positions of major constellations for the sky map"""
        self.constellations = {
            'Orion': {'ra': 83.8, 'dec': -5.4, 'size': 30},
            'Ursa Major': {'ra': 165.0, 'dec': 55.0, 'size': 25},
            'Cassiopeia': {'ra': 15.0, 'dec': 60.0, 'size': 20},
            'Andromeda': {'ra': 10.7, 'dec': 41.3, 'size': 20},
            'Cygnus': {'ra': 305.0, 'dec': 40.0, 'size': 25},
            'Sagittarius': {'ra': 283.0, 'dec': -25.0, 'size': 30},
            'Taurus': {'ra': 68.0, 'dec': 15.0, 'size': 25},
            'Leo': {'ra': 152.0, 'dec': 15.0, 'size': 25},
            'Aquila': {'ra': 297.0, 'dec': 10.0, 'size': 20},
            'Centaurus': {'ra': 192.0, 'dec': -43.0, 'size': 30},
        }
        
        # Famous objects
        self.famous_objects = {
            'Horsehead': {'ra': 85.2479, 'dec': -2.4583, 'label': 'Horsehead'},
            'M31': {'ra': 10.6847, 'dec': 41.2687, 'label': 'M31 (Andromeda)'},
            'M42': {'ra': 83.8221, 'dec': -5.3911, 'label': 'M42 (Orion)'},
            'M51': {'ra': 202.4696, 'dec': 47.1952, 'label': 'M51 (Whirlpool)'},
            'M57': {'ra': 283.3963, 'dec': 33.0295, 'label': 'M57 (Ring)'},
            'Crab': {'ra': 83.6333, 'dec': 22.0145, 'label': 'M1 (Crab)'},
        }
        
    def setup_ui(self):
        """Create the user interface"""
        
        # Main container with three panels
        main_frame = ttk.Frame(self.root, padding="10")
        main_frame.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
        
        self.root.columnconfigure(0, weight=1)
        self.root.rowconfigure(0, weight=1)
        main_frame.columnconfigure(1, weight=2)
        main_frame.columnconfigure(2, weight=3)
        main_frame.rowconfigure(0, weight=1)
        
        # LEFT PANEL - Sky Map
        sky_frame = ttk.LabelFrame(main_frame, text="🌌 Interactive Sky Map - Click Anywhere!", padding="10")
        sky_frame.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S), padx=(0, 5))
        sky_frame.columnconfigure(0, weight=1)
        sky_frame.rowconfigure(0, weight=1)
        
        # Create sky map figure
        self.sky_fig = Figure(figsize=(6, 6), dpi=100)
        self.sky_ax = self.sky_fig.add_subplot(111)
        self.sky_canvas = FigureCanvasTkAgg(self.sky_fig, master=sky_frame)
        self.sky_canvas.get_tk_widget().grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
        
        # Connect click event
        self.sky_canvas.mpl_connect('button_press_event', self.on_sky_click)
        
        # Sky map info
        sky_info = ttk.Label(sky_frame, text="Click on the sky map to select coordinates\n" +
                            "Constellations and famous objects are marked\n" +
                            "RA: 0-360°, Dec: -90° to +90°", 
                            justify=tk.CENTER)
        sky_info.grid(row=1, column=0, pady=5)
        
        # MIDDLE PANEL - Controls
        control_frame = ttk.LabelFrame(main_frame, text="⚙️ Controls", padding="10")
        control_frame.grid(row=0, column=1, sticky=(tk.W, tk.E, tk.N, tk.S), padx=5)
        
        # Selected coordinates display
        ttk.Label(control_frame, text="Selected Coordinates:", 
                 font=('TkDefaultFont', 10, 'bold')).grid(row=0, column=0, sticky=tk.W, pady=(0, 10))
        
        coord_display = ttk.Frame(control_frame)
        coord_display.grid(row=1, column=0, sticky=(tk.W, tk.E), pady=5)
        
        ttk.Label(coord_display, text="RA:").grid(row=0, column=0, sticky=tk.W)
        self.ra_display = ttk.Label(coord_display, text="---.----°", 
                                    foreground="blue", font=('TkDefaultFont', 11, 'bold'))
        self.ra_display.grid(row=0, column=1, sticky=tk.W, padx=5)
        
        ttk.Label(coord_display, text="Dec:").grid(row=1, column=0, sticky=tk.W)
        self.dec_display = ttk.Label(coord_display, text="---.----°", 
                                     foreground="blue", font=('TkDefaultFont', 11, 'bold'))
        self.dec_display.grid(row=1, column=1, sticky=tk.W, padx=5)
        
        # Manual entry option
        ttk.Separator(control_frame, orient='horizontal').grid(row=2, column=0, sticky=(tk.W, tk.E), pady=15)
        
        ttk.Label(control_frame, text="Or Enter Manually:", 
                 font=('TkDefaultFont', 10, 'bold')).grid(row=3, column=0, sticky=tk.W, pady=5)
        
        ttk.Label(control_frame, text="RA (0-360°):").grid(row=4, column=0, sticky=tk.W)
        self.ra_entry = ttk.Entry(control_frame, width=20)
        self.ra_entry.grid(row=5, column=0, sticky=(tk.W, tk.E), pady=2)
        
        ttk.Label(control_frame, text="Dec (-90 to +90°):").grid(row=6, column=0, sticky=tk.W)
        self.dec_entry = ttk.Entry(control_frame, width=20)
        self.dec_entry.grid(row=7, column=0, sticky=(tk.W, tk.E), pady=2)
        
        ttk.Button(control_frame, text="Set Coordinates", 
                  command=self.set_manual_coords).grid(row=8, column=0, sticky=(tk.W, tk.E), pady=5)
        
        # Observation parameters
        ttk.Separator(control_frame, orient='horizontal').grid(row=9, column=0, sticky=(tk.W, tk.E), pady=15)
        
        ttk.Label(control_frame, text="Observation Settings:", 
                 font=('TkDefaultFont', 10, 'bold')).grid(row=10, column=0, sticky=tk.W, pady=5)
        
        ttk.Label(control_frame, text="Field of View:").grid(row=11, column=0, sticky=tk.W)
        self.fov_var = tk.StringVar(value="0.5")
        fov_frame = ttk.Frame(control_frame)
        fov_frame.grid(row=12, column=0, sticky=(tk.W, tk.E))
        ttk.Entry(fov_frame, textvariable=self.fov_var, width=8).pack(side=tk.LEFT)
        ttk.Label(fov_frame, text="° (degrees)").pack(side=tk.LEFT, padx=5)
        
        # Quick FOV buttons
        fov_quick = ttk.Frame(control_frame)
        fov_quick.grid(row=13, column=0, sticky=(tk.W, tk.E), pady=2)
        ttk.Button(fov_quick, text="0.1°", width=5,
                  command=lambda: self.fov_var.set("0.1")).pack(side=tk.LEFT, padx=2)
        ttk.Button(fov_quick, text="0.5°", width=5,
                  command=lambda: self.fov_var.set("0.5")).pack(side=tk.LEFT, padx=2)
        ttk.Button(fov_quick, text="1.0°", width=5,
                  command=lambda: self.fov_var.set("1.0")).pack(side=tk.LEFT, padx=2)
        ttk.Button(fov_quick, text="2.0°", width=5,
                  command=lambda: self.fov_var.set("2.0")).pack(side=tk.LEFT, padx=2)
        
        ttk.Label(control_frame, text="Survey:").grid(row=14, column=0, sticky=tk.W, pady=(10, 0))
        self.survey_var = tk.StringVar(value="DSS2 Red")
        surveys = ['DSS2 Red', 'DSS2 Blue', 'DSS2 IR', '2MASS-J', '2MASS-H', '2MASS-K', 'WISE 3.4']
        ttk.Combobox(control_frame, textvariable=self.survey_var, 
                    values=surveys, width=18, state='readonly').grid(row=15, column=0, sticky=(tk.W, tk.E), pady=2)
        
        ttk.Label(control_frame, text="Resolution:").grid(row=16, column=0, sticky=tk.W)
        self.resolution_var = tk.StringVar(value="1000")
        ttk.Combobox(control_frame, textvariable=self.resolution_var, 
                    values=['500', '800', '1000', '1500', '2000'], 
                    width=18, state='readonly').grid(row=17, column=0, sticky=(tk.W, tk.E), pady=2)
        
        # Download button
        ttk.Button(control_frame, text="📡 Download & View", 
                  command=self.download_and_view,
                  style='Accent.TButton').grid(row=18, column=0, sticky=(tk.W, tk.E), pady=20)
        
        # Display options
        ttk.Label(control_frame, text="Display:").grid(row=19, column=0, sticky=tk.W)
        self.cmap_var = tk.StringVar(value="gray")
        ttk.Combobox(control_frame, textvariable=self.cmap_var, 
                    values=['gray', 'viridis', 'hot', 'cool', 'rainbow', 'plasma', 'twilight'],
                    width=18, state='readonly').grid(row=20, column=0, sticky=(tk.W, tk.E), pady=2)
        
        ttk.Button(control_frame, text="Update View", 
                  command=self.update_display).grid(row=21, column=0, sticky=(tk.W, tk.E), pady=2)
        
        # Save buttons
        ttk.Button(control_frame, text="💾 Save FITS", 
                  command=self.save_fits).grid(row=22, column=0, sticky=(tk.W, tk.E), pady=2)
        ttk.Button(control_frame, text="💾 Save PNG", 
                  command=self.save_png).grid(row=23, column=0, sticky=(tk.W, tk.E), pady=2)
        
        # Status
        self.status_label = ttk.Label(control_frame, text="Ready - Click on sky map", foreground="green")
        self.status_label.grid(row=24, column=0, sticky=tk.W, pady=10)
        
        # RIGHT PANEL - Image Display
        display_frame = ttk.LabelFrame(main_frame, text="🔭 Telescope View", padding="10")
        display_frame.grid(row=0, column=2, sticky=(tk.W, tk.E, tk.N, tk.S), padx=(5, 0))
        display_frame.columnconfigure(0, weight=1)
        display_frame.rowconfigure(1, weight=1)
        
        # Info text
        self.info_text = tk.Text(display_frame, height=3, width=60, wrap=tk.WORD)
        self.info_text.grid(row=0, column=0, sticky=(tk.W, tk.E), pady=(0, 10))
        self.info_text.insert('1.0', "Click anywhere on the sky map to select coordinates.\n" +
                             "Then click 'Download & View' to observe that region.")
        self.info_text.config(state='disabled')
        
        # Create matplotlib figure for telescope view
        self.view_fig = Figure(figsize=(8, 7), dpi=100)
        self.view_ax = self.view_fig.add_subplot(111)
        self.view_canvas = FigureCanvasTkAgg(self.view_fig, master=display_frame)
        self.view_canvas.get_tk_widget().grid(row=1, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
        
        # Add matplotlib toolbar in a container frame (to avoid pack/grid conflict)
        toolbar_frame = ttk.Frame(display_frame)
        toolbar_frame.grid(row=2, column=0, sticky=(tk.W, tk.E))
        toolbar = NavigationToolbar2Tk(self.view_canvas, toolbar_frame)
        toolbar.update()
        
        # Initial display
        self.show_welcome_view()
        
    def draw_sky_map(self):
        """Draw the interactive all-sky map"""
        self.sky_ax.clear()
        
        # Set up Aitoff projection (all-sky map)
        self.sky_ax.set_xlim(360, 0)  # RA increases to the left (astronomical convention)
        self.sky_ax.set_ylim(-90, 90)
        
        # Grid
        self.sky_ax.grid(True, alpha=0.3, linestyle='--')
        self.sky_ax.set_xlabel('Right Ascension (degrees)', fontsize=10)
        self.sky_ax.set_ylabel('Declination (degrees)', fontsize=10)
        self.sky_ax.set_title('All-Sky Map (Click to Select)', fontsize=12, fontweight='bold')
        
        # Draw equator
        ra_line = np.linspace(0, 360, 100)
        self.sky_ax.plot(ra_line, np.zeros_like(ra_line), 'b-', alpha=0.3, linewidth=2, label='Celestial Equator')
        
        # Draw ecliptic (approximate)
        ecliptic_dec = 23.4 * np.sin(np.radians(ra_line))
        self.sky_ax.plot(ra_line, ecliptic_dec, 'orange', alpha=0.3, linewidth=1, linestyle='--', label='Ecliptic')
        
        # Draw Milky Way plane (very approximate)
        # Galactic center is around RA~266°, Dec~-29°
        mw_ra = np.linspace(0, 360, 100)
        mw_dec = -29 + 60 * np.sin(np.radians(mw_ra - 266))
        self.sky_ax.plot(mw_ra, mw_dec, 'purple', alpha=0.2, linewidth=3, label='Milky Way')
        
        # Plot constellations
        for name, data in self.constellations.items():
            circle = plt.Circle((data['ra'], data['dec']), data['size']/2, 
                               color='cyan', alpha=0.1, zorder=1)
            self.sky_ax.add_patch(circle)
            self.sky_ax.text(data['ra'], data['dec'], name, 
                           ha='center', va='center', fontsize=8, 
                           alpha=0.6, color='cyan', weight='bold')
        
        # Plot famous objects
        for name, data in self.famous_objects.items():
            self.sky_ax.plot(data['ra'], data['dec'], 'r*', markersize=12, zorder=3)
            self.sky_ax.text(data['ra'], data['dec'] + 5, data['label'], 
                           ha='center', va='bottom', fontsize=7, 
                           color='red', weight='bold')
        
        # Plot selected point if exists
        if self.selected_ra is not None and self.selected_dec is not None:
            self.sky_ax.plot(self.selected_ra, self.selected_dec, 
                           'go', markersize=15, markeredgewidth=2, 
                           markeredgecolor='yellow', zorder=5)
            self.sky_ax.text(self.selected_ra, self.selected_dec - 8, 
                           'SELECTED', ha='center', va='top', 
                           fontsize=9, color='yellow', weight='bold',
                           bbox=dict(boxstyle='round', facecolor='green', alpha=0.7))
        
        self.sky_ax.legend(loc='upper right', fontsize=8)
        self.sky_ax.invert_xaxis()  # RA increases to the left
        
        self.sky_canvas.draw()
    
    def on_sky_click(self, event):
        """Handle click on sky map"""
        if event.inaxes == self.sky_ax and event.xdata and event.ydata:
            self.selected_ra = event.xdata
            self.selected_dec = event.ydata
            
            # Update display
            self.ra_display.config(text=f"{self.selected_ra:.4f}°")
            self.dec_display.config(text=f"{self.selected_dec:.4f}°")
            
            # Update entry fields
            self.ra_entry.delete(0, tk.END)
            self.ra_entry.insert(0, f"{self.selected_ra:.4f}")
            self.dec_entry.delete(0, tk.END)
            self.dec_entry.insert(0, f"{self.selected_dec:.4f}")
            
            # Redraw sky map with selection marker
            self.draw_sky_map()
            
            # Convert to HMS/DMS for info
            coord = SkyCoord(ra=self.selected_ra*u.degree, dec=self.selected_dec*u.degree)
            ra_hms = coord.ra.to_string(unit=u.hour, sep=':', precision=1)
            dec_dms = coord.dec.to_string(unit=u.degree, sep=':', precision=1)
            
            info = f"Selected: RA={self.selected_ra:.4f}° ({ra_hms}), Dec={self.selected_dec:.4f}° ({dec_dms})\n"
            info += "Click 'Download & View' to observe this region."
            self.update_info(info)
            
            self.status_label.config(text="Coordinates selected - Ready to download", foreground="blue")
    
    def set_manual_coords(self):
        """Set coordinates from manual entry"""
        try:
            ra = float(self.ra_entry.get())
            dec = float(self.dec_entry.get())
            
            if not (0 <= ra <= 360):
                messagebox.showerror("Invalid RA", "Right Ascension must be between 0 and 360 degrees")
                return
            if not (-90 <= dec <= 90):
                messagebox.showerror("Invalid Dec", "Declination must be between -90 and +90 degrees")
                return
            
            self.selected_ra = ra
            self.selected_dec = dec
            
            self.ra_display.config(text=f"{self.selected_ra:.4f}°")
            self.dec_display.config(text=f"{self.selected_dec:.4f}°")
            
            self.draw_sky_map()
            
            self.status_label.config(text="Coordinates set - Ready to download", foreground="blue")
            
        except ValueError:
            messagebox.showerror("Invalid Input", "Please enter valid numerical coordinates")
    
    def download_and_view(self):
        """Download image for selected coordinates"""
        if self.selected_ra is None or self.selected_dec is None:
            messagebox.showwarning("No Selection", "Please select coordinates on the sky map first")
            return
        
        try:
            fov = float(self.fov_var.get())
            survey = self.survey_var.get()
            resolution = int(self.resolution_var.get())
            
            self.status_label.config(text="Downloading from SkyView...", foreground="orange")
            self.root.update()
            
            # Run in thread
            thread = threading.Thread(target=self._download_thread, 
                                     args=(self.selected_ra, self.selected_dec, fov, survey, resolution))
            thread.daemon = True
            thread.start()
            
        except ValueError as e:
            messagebox.showerror("Invalid Input", f"Please check your values:\n{e}")
    
    def _download_thread(self, ra, dec, fov, survey, resolution):
        """Background download"""
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
                
                with fits.open(BytesIO(self.current_fits_data)) as hdul:
                    self.current_image_data = hdul[0].data
                    header = hdul[0].header
                
                self.root.after(0, self._display_image, header, survey, ra, dec, fov)
                self.root.after(0, lambda: self.status_label.config(text="✓ Download complete!", foreground="green"))
            else:
                self.root.after(0, lambda: messagebox.showerror("Download Failed", 
                               f"Failed to retrieve image. Status: {response.status_code}"))
                self.root.after(0, lambda: self.status_label.config(text="✗ Download failed", foreground="red"))
                
        except Exception as e:
            self.root.after(0, lambda: messagebox.showerror("Error", str(e)))
            self.root.after(0, lambda: self.status_label.config(text="✗ Error occurred", foreground="red"))
    
    def _display_image(self, header, survey, ra, dec, fov):
        """Display downloaded image"""
        # Clear the figure completely to remove old colorbars
        self.view_fig.clear()
        self.view_ax = self.view_fig.add_subplot(111)
        
        vmin = np.nanpercentile(self.current_image_data, 1)
        vmax = np.nanpercentile(self.current_image_data, 99)
        
        im = self.view_ax.imshow(self.current_image_data, 
                                cmap=self.cmap_var.get(), 
                                origin='lower',
                                vmin=vmin, 
                                vmax=vmax,
                                interpolation='nearest')
        
        self.view_fig.colorbar(im, ax=self.view_ax, label='Intensity (ADU)', fraction=0.046)
        
        # Convert to HMS/DMS
        coord = SkyCoord(ra=ra*u.degree, dec=dec*u.degree)
        ra_hms = coord.ra.to_string(unit=u.hour, sep=':', precision=1)
        dec_dms = coord.dec.to_string(unit=u.degree, sep=':', precision=1)
        
        self.view_ax.set_title(f'{survey} - FOV: {fov}°\nRA: {ra_hms} ({ra:.4f}°), Dec: {dec_dms} ({dec:.4f}°)', 
                              fontsize=11)
        self.view_ax.set_xlabel('X (pixels)')
        self.view_ax.set_ylabel('Y (pixels)')
        
        self.view_canvas.draw()
        
        # Update info
        info = f"Survey: {survey} | FOV: {fov}° | Resolution: {self.current_image_data.shape}\n"
        info += f"Position: RA={ra:.4f}° ({ra_hms}), Dec={dec:.4f}° ({dec_dms})\n"
        info += f"Intensity: [{np.nanmin(self.current_image_data):.2e}, {np.nanmax(self.current_image_data):.2e}] ADU"
        self.update_info(info)
    
    def update_display(self):
        """Update colormap"""
        if self.current_image_data is not None:
            # Clear figure to remove old colorbars
            self.view_fig.clear()
            self.view_ax = self.view_fig.add_subplot(111)
            
            ra = self.selected_ra
            dec = self.selected_dec
            fov = float(self.fov_var.get())
            survey = self.survey_var.get()
            
            vmin = np.nanpercentile(self.current_image_data, 1)
            vmax = np.nanpercentile(self.current_image_data, 99)
            
            im = self.view_ax.imshow(self.current_image_data, 
                                   cmap=self.cmap_var.get(), 
                                   origin='lower',
                                   vmin=vmin, 
                                   vmax=vmax,
                                   interpolation='nearest')
            
            self.view_fig.colorbar(im, ax=self.view_ax, label='Intensity (ADU)', fraction=0.046)
            
            # Convert to HMS/DMS
            coord = SkyCoord(ra=ra*u.degree, dec=dec*u.degree)
            ra_hms = coord.ra.to_string(unit=u.hour, sep=':', precision=1)
            dec_dms = coord.dec.to_string(unit=u.degree, sep=':', precision=1)
            
            self.view_ax.set_title(f'{survey} - FOV: {fov}°\nRA: {ra_hms} ({ra:.4f}°), Dec: {dec_dms} ({dec:.4f}°)', 
                                  fontsize=11)
            self.view_ax.set_xlabel('X (pixels)')
            self.view_ax.set_ylabel('Y (pixels)')
            
            self.view_canvas.draw()
        else:
            messagebox.showinfo("No Image", "Download an image first")
    
    def save_fits(self):
        """Save FITS file"""
        if self.current_fits_data is None:
            messagebox.showinfo("No Image", "Download an image first")
            return
        
        filename = filedialog.asksaveasfilename(
            defaultextension=".fits",
            filetypes=[("FITS files", "*.fits"), ("All files", "*.*")],
            initialfile=f"sky_RA{self.selected_ra:.2f}_Dec{self.selected_dec:.2f}.fits"
        )
        
        if filename:
            with open(filename, 'wb') as f:
                f.write(self.current_fits_data)
            messagebox.showinfo("Saved", f"FITS saved:\n{filename}")
    
    def save_png(self):
        """Save PNG"""
        if self.current_image_data is None:
            messagebox.showinfo("No Image", "Download an image first")
            return
        
        filename = filedialog.asksaveasfilename(
            defaultextension=".png",
            filetypes=[("PNG files", "*.png"), ("All files", "*.*")],
            initialfile=f"sky_RA{self.selected_ra:.2f}_Dec{self.selected_dec:.2f}.png"
        )
        
        if filename:
            self.view_fig.savefig(filename, dpi=150, bbox_inches='tight')
            messagebox.showinfo("Saved", f"PNG saved:\n{filename}")
    
    def show_welcome_view(self):
        """Initial view display"""
        self.view_ax.clear()
        self.view_ax.text(0.5, 0.5, 
                         '🔭 Telescope View\n\n' +
                         'Click on the sky map to select\n' +
                         'any region of the sky,\n' +
                         'then download to view it here',
                         ha='center', va='center', fontsize=14, 
                         transform=self.view_ax.transAxes)
        self.view_ax.axis('off')
        self.view_canvas.draw()
    
    def update_info(self, text):
        """Update info panel"""
        self.info_text.config(state='normal')
        self.info_text.delete('1.0', tk.END)
        self.info_text.insert('1.0', text)
        self.info_text.config(state='disabled')

def main():
    root = tk.Tk()
    app = InteractiveSkyViewer(root)
    root.mainloop()

if __name__ == "__main__":
    main()
