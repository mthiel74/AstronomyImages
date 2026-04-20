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

NASA SkyView aggregates data from many surveys. **Licensing depends on who produced the survey, not where you downloaded it from.** The two regimes are:

#### NASA-operated surveys (2MASS, WISE, GALEX, ROSAT, NVSS, …)
- **License**: Public Domain (NASA work, 17 U.S.C. § 105; or joint NASA/NSF mission releases)
- **Policy**: https://skyview.gsfc.nasa.gov/about.html
- **What you CAN do**:
  - ✅ Use for research, education, or personal projects
  - ✅ Post images on GitHub, LinkedIn, Twitter, or any social media
  - ✅ Include in presentations, papers, and publications
  - ✅ Use commercially without restrictions
- **Attribution (recommended, not required)** for 2MASS: *"This publication makes use of data products from the Two Micron All Sky Survey, which is a joint project of the University of Massachusetts and the Infrared Processing and Analysis Center / California Institute of Technology, funded by the National Aeronautics and Space Administration and the National Science Foundation."*
- **Attribution (recommended, not required)** for WISE: *"This publication makes use of data products from the Wide-field Infrared Survey Explorer, a joint project of UCLA and JPL/Caltech, funded by NASA."*

#### DSS / DSS2 (Digitized Sky Survey, STScI-distributed)
- **License**: **Not public domain.** Free for research and educational use; commercial use requires contacting STScI.
- **Copyright**: held by the Anglo-Australian Observatory Board, Caltech (Palomar plates), AURA, and UK SERC/PPARC + the Anglo-Australian Telescope Board, depending on the plate. Distributed by STScI under research/education terms.
- **Policy**: https://archive.stsci.edu/dss/acknowledging.html (+ https://archive.stsci.edu/dss/copyright.html)
- **What you CAN do**:
  - ✅ Use for academic research and education
  - ✅ Post on websites, blogs, and social media **with the attribution text below**
  - ✅ Include in publications with the attribution text below
  - ⚠️ Commercial use: contact STScI
- **Attribution Required (verbatim)**: *"The Digitized Sky Surveys were produced at the Space Telescope Science Institute under U.S. Government grant NAG W-2166. The images of these surveys are based on photographic data obtained using the Oschin Schmidt Telescope on Palomar Mountain and the UK Schmidt Telescope."*

### Summary
- NASA-operated surveys (2MASS, WISE, GALEX, …): **public domain** — use freely, attribution optional.
- DSS / DSS2: **copyrighted, research/education use with required STScI attribution**; commercial use requires contacting STScI.
