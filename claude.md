# AstronomyImages

Interactive astronomy image viewer and analysis toolkit for accessing real telescope data and virtual observatory archives.

## Project Overview

This project provides Python-based tools for downloading, visualizing, and analyzing astronomical images from public sky surveys and virtual observatories. It includes both command-line scripts and interactive GUI applications for exploring the night sky.

## Key Features

- **Interactive Sky Map Viewer**: Click-based interface to explore any region of the sky with real telescope data
- **Multiple Survey Support**: Access to DSS2 (Red/Blue/IR), 2MASS (J/H/K bands), WISE, and other sky surveys
- **Real-Time Downloads**: Fetch FITS astronomical images via NASA SkyView API
- **Advanced Visualization**: Multiple color maps, intensity scaling, and analysis tools
- **FITS File Support**: Full support for astronomical FITS format with WCS (World Coordinate System)
- **Simulation Capabilities**: Generate simulated observations for educational purposes (Horsehead Nebula demo)

## Main Scripts

### Interactive Applications
- **`interactive_sky_map.py`**: Full-featured GUI with clickable all-sky map, constellation overlays, and telescope view
- **`astronomy_viewer_gui.py`**: Simplified viewer with object catalog and manual coordinate entry

### Demonstration & Examples
- **`horsehead_complete_demo.py`**: Complete workflow demonstration with simulated Horsehead Nebula observation
  - Physics-based simulation of dark nebulae
  - Optical depth calculations
  - Comprehensive 9-panel analysis plots
  - API documentation and examples

### Utility Scripts
- **`horsehead_skyview.py`** / **`horsehead_skyview_fixed.py`**: Direct SkyView API integration examples
- **`horsehead_alternative.py`**: Alternative implementation approaches

## Data Access APIs

### NASA SkyView Virtual Observatory
- **URL**: https://skyview.gsfc.nasa.gov
- **Authentication**: None required
- **Surveys**: DSS2, 2MASS, WISE, SDSS, and 30+ others
- **Format**: FITS images with WCS headers

### Las Cumbres Observatory (Educational Access)
- **URL**: https://observe.lco.global
- **Authentication**: Free API token for education
- **Capability**: Schedule real telescope observations worldwide
- **Details**: See `api_examples.json` for request format

## File Structure

```
AstronomyImages/
├── interactive_sky_map.py          # Main interactive viewer
├── astronomy_viewer_gui.py         # Simplified GUI
├── horsehead_complete_demo.py      # Full demo with simulation
├── horsehead_skyview*.py           # API integration examples
├── api_examples.json               # API documentation
├── *.fits                          # Downloaded astronomical images
├── *.png                           # Visualization outputs
└── horsehead_analysis.json         # Analysis results
```

## Dependencies

```python
# Core astronomy libraries
astropy              # FITS files, WCS, coordinates
numpy                # Numerical operations
matplotlib           # Visualization

# GUI (for interactive viewers)
tkinter              # Built-in Python GUI

# Network
requests             # API calls
```

## Usage Examples

### Interactive Sky Map
```bash
python interactive_sky_map.py
# Click anywhere on the sky map to select coordinates
# Adjust FOV, survey, and resolution
# Download and view in real-time
```

### Simple Viewer
```bash
python astronomy_viewer_gui.py
# Choose from catalog of famous objects
# Or enter manual coordinates
# Download FITS and save locally
```

### Run Horsehead Demo
```bash
python horsehead_complete_demo.py
# Generates simulated observation
# Creates comprehensive analysis plots
# Outputs FITS, PNG, and JSON files
```

## Technical Notes

- All coordinates use J2000 equatorial system (RA/Dec in degrees)
- FITS images include WCS headers for coordinate mapping
- Default pixel scale: 1 arcsec/pixel (configurable)
- Image scaling: Linear or percentile-based (1-99%)
- No API keys required for SkyView access

## Educational Applications

- Learn about astronomical coordinate systems
- Explore real telescope data interactively
- Understand FITS file format and WCS
- Practice image analysis techniques
- Visualize dark nebulae and emission regions

## Security & Privacy

- **No authentication credentials stored in code**
- All API tokens are placeholder examples only ("YOUR_API_TOKEN")
- Safe for public repository sharing
- Network requests use standard HTTPS

## Future Enhancements

- Add support for more sky surveys (Pan-STARRS, etc.)
- Implement advanced image processing (stacking, calibration)
- Integration with SIMBAD/NED for object identification
- Export to DS9 region files
- Multi-wavelength comparison views

## License & Data Usage Rights

### This Software (Python Code)
**MIT License** - Free to use, modify, and distribute for any purpose, including commercial applications.

Copyright (c) 2025 - This astronomy viewer toolkit

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software.

### Astronomical Survey Data Licenses

All astronomical survey data accessed through this tool comes from public archives with permissive usage terms:

#### NASA SkyView Virtual Observatory
- **License**: Public Domain (U.S. Government Work, 17 U.S.C. § 105)
- **Official Policy**: https://skyview.gsfc.nasa.gov/about.html
- **Surveys Covered**: DSS2 Red/Blue/IR, WISE 3.4/4.6/12/22, GALEX, RASS, NVSS, and all NASA-operated surveys
- **Usage Rights**:
  - ✅ Unrestricted use for any purpose (research, education, commercial)
  - ✅ No permission required for publication, redistribution, or derivative works
  - ✅ May be posted on websites, social media (GitHub, LinkedIn, Twitter, etc.)
  - ✅ May be included in papers, presentations, books, or courses
  - ✅ May be sold or used in commercial products
- **Attribution**: Not legally required, but professionally courteous
  - Suggested: "Image courtesy of NASA/SkyView"

#### 2MASS (Two Micron All Sky Survey)
- **License**: Public Domain
- **Official Policy**: https://www.ipac.caltech.edu/2mass/overview/about2mass.html
- **Data Rights Policy**: https://www.ipac.caltech.edu/2mass/releases/allsky/doc/sec1_6b.html
- **Surveys**: 2MASS-J (1.25 µm), 2MASS-H (1.65 µm), 2MASS-K (2.17 µm)
- **Usage Rights**: Identical to NASA data (public domain, unrestricted use)
- **Recommended Attribution**:
  - "This publication makes use of data products from the Two Micron All Sky Survey, which is a joint project of the University of Massachusetts and the Infrared Processing and Analysis Center/California Institute of Technology, funded by NASA and NSF."

#### DSS (Digitized Sky Survey) via STScI
- **License**: Free for research and educational purposes
- **Official Policy**: https://archive.stsci.edu/dss/acknowledging.html
- **Surveys**: DSS1, DSS2 (via SkyView access)
- **Usage Rights**:
  - ✅ Free use for academic research and education
  - ✅ May be posted on personal/institutional websites with attribution
  - ✅ May be included in publications with proper citation
  - ✅ May be shared on social media (GitHub, LinkedIn, portfolio sites)
- **Required Attribution**:
  - "The Digitized Sky Survey was produced at the Space Telescope Science Institute under U.S. Government grant NAG W-2166. The images of these surveys are based on photographic data obtained using the Oschin Schmidt Telescope on Palomar Mountain and the UK Schmidt Telescope."

#### WISE (Wide-field Infrared Survey Explorer)
- **License**: Public Domain (NASA data)
- **Official Policy**: https://wise2.ipac.caltech.edu/docs/release/allsky/
- **Usage Rights**: Unrestricted (public domain)
- **Recommended Attribution**: "Data from WISE (Wide-field Infrared Survey Explorer), a joint project of UCLA and JPL/Caltech, funded by NASA"

### API Access Terms

#### NASA SkyView API
- **Endpoint**: https://skyview.gsfc.nasa.gov/cgi-bin/images
- **Rate Limits**: "Reasonable use" (no hard limit published; avoid automated bulk downloads)
- **Terms**: https://skyview.gsfc.nasa.gov/about.html
- **Contact**: skyview@skyview.gsfc.nasa.gov for bulk access arrangements

#### Las Cumbres Observatory API
- **Educational Access**: Free API tokens available at https://observe.lco.global/
- **Terms**: https://lco.global/observatory/data/
- **Usage**: Educational/research observations are free; commercial use requires separate agreement
- **Data Rights**: You own copyright to observations you schedule

### Can I Post Images on Social Media?

**YES!** Here's what you can do:

| Content Type | GitHub | LinkedIn | Twitter | Personal Website | Commercial Use |
|--------------|--------|----------|---------|------------------|----------------|
| NASA/SkyView images | ✅ Yes | ✅ Yes | ✅ Yes | ✅ Yes | ✅ Yes |
| 2MASS images | ✅ Yes | ✅ Yes | ✅ Yes | ✅ Yes | ✅ Yes |
| DSS images (with attribution) | ✅ Yes | ✅ Yes | ✅ Yes | ✅ Yes | ⚠️ Contact STScI |

### Best Practices for Attribution

When sharing astronomical images publicly:

**Minimal (Social Media)**:
```
Image: NASA SkyView / Horsehead Nebula (DSS2)
```

**Standard (Blog/Portfolio)**:
```
Astronomical image retrieved from NASA SkyView Virtual Observatory
Survey: Digitized Sky Survey 2 (DSS2), Red filter
Object: Horsehead Nebula (Barnard 33)
```

**Formal (Publications)**:
```
This research has made use of the NASA/GSFC SkyView Virtual Observatory
(https://skyview.gsfc.nasa.gov). DSS2 images were produced at STScI under
U.S. Government grant NAG W-2166.
```

### Summary

**For users of this software**:
- 🟢 All NASA and 2MASS data is **public domain** - use freely with no restrictions
- 🟢 DSS data is free for education/research - just add attribution
- 🟢 The Python code is **MIT licensed** - modify and redistribute freely
- 🟢 **Safe to post on GitHub, LinkedIn, and all social media**

**No confidential data or API keys are included in this repository.**
