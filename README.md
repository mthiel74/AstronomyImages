# AstronomyImages 🔭

Interactive astronomy image viewer and analysis toolkit for accessing real telescope data from virtual observatories.

![AstronomyImages Interface](AstroInterface.png)

## Quick Start

```bash
# Install dependencies
pip install astropy numpy matplotlib requests

# Run the interactive sky map
python interactive_sky_map.py

# Or run the simplified viewer
python astronomy_viewer_gui.py
```

## Features

- 🌌 **Interactive Sky Map**: Click anywhere on an all-sky map to explore that region
- 🔭 **Real Telescope Data**: Download FITS images from NASA SkyView and other archives
- 📊 **Multiple Surveys**: DSS2, 2MASS, WISE, and 30+ sky surveys
- 🎨 **Advanced Visualization**: Multiple colormaps, scaling options, and analysis tools
- 💾 **Export**: Save FITS files and publication-quality PNG images

## Main Applications

### Interactive Sky Map (`interactive_sky_map.py`)
Full-featured GUI with:
- Clickable all-sky map showing constellations and famous objects
- Real-time coordinate selection
- Adjustable field of view and resolution
- Multiple survey options

### Astronomy Viewer (`astronomy_viewer_gui.py`)
Simplified interface with:
- Catalog of popular astronomical objects
- Manual coordinate entry
- Quick download and visualization

### Horsehead Nebula Demo (`horsehead_complete_demo.py`)
Educational demonstration featuring:
- Physics-based simulation of dark nebulae
- Comprehensive 9-panel analysis
- API usage examples
- FITS file creation with WCS headers

## No API Keys Required

All public data access is free and requires no authentication. The SkyView API is open for educational and research use.

## Documentation

See `claude.md` for complete technical documentation, API details, and usage examples.

## License & Data Usage

### This Software (Python Code)
**MIT License** - Free to use, modify, and distribute for any purpose, including commercial applications.

### Astronomical Data Usage Rights

All astronomical survey data accessed through this tool is **publicly available** and can be used freely with proper attribution:

#### NASA SkyView Data (DSS2, WISE, etc.)
- **License**: Public Domain (U.S. Government Work)
- **Policy**: https://skyview.gsfc.nasa.gov/about.html
- **What you CAN do**:
  - ✅ Use for research, education, or personal projects
  - ✅ Post images on GitHub, LinkedIn, Twitter, or any social media
  - ✅ Include in presentations, papers, and publications
  - ✅ Use commercially without restrictions
- **Attribution**: Recommended but not required. Example: "Image courtesy of SkyView/NASA"

#### 2MASS Survey Data
- **License**: Public Domain
- **Policy**: https://www.ipac.caltech.edu/2mass/overview/about2mass.html
- **What you CAN do**: Same as NASA data (public domain)
- **Attribution**: "This publication makes use of data products from the Two Micron All Sky Survey (2MASS)"

#### DSS (Digitized Sky Survey)
- **License**: Free for research and educational use
- **Policy**: https://archive.stsci.edu/dss/acknowledging.html
- **What you CAN do**:
  - ✅ Use for academic research and education
  - ✅ Post on websites and social media with attribution
  - ✅ Include in publications with citation
- **Attribution Required**: "The Digitized Sky Survey was produced at the Space Telescope Science Institute"

### Summary
All astronomical data from SkyView and 2MASS is **public domain** - use it however you want! DSS data requires simple attribution.
