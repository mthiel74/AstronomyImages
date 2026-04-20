# Wolfram Language Port

This directory contains Wolfram Language translations of the Python scripts in the parent `AstronomyImages/` project. Each `.wl` file is a standalone package wrapped in `BeginPackage` / `EndPackage` so it can be loaded with `Get` or run directly from the command line via `wolframscript`.

## File mapping

| Python source | Wolfram port | What it does |
|---|---|---|
| `api_examples.json` | `APIExamples.wl` | Stores SkyView + LCO request templates as an `Association` and exports them back to JSON. |
| `horsehead_skyview.py` | `HorseheadSkyView.wl` | Queries NASA SkyView for DSS2 / 2MASS cutouts of the Horsehead Nebula and saves FITS + PNG. |
| `horsehead_skyview_fixed.py` | `HorseheadSkyViewFixed.wl` | Same, but writes to the current directory. |
| `horsehead_alternative.py` | `HorseheadAlternative.wl` | SDSS SkyServer JPEG + STScI DSS2 JPEG fallbacks with histogram analysis. |
| `horsehead_complete_demo.py` | `HorseheadCompleteDemo.wl` | Physics-based simulation, FITS writer with WCS-style header, 9-panel analysis. |
| `astronomy_viewer_gui.py` | `AstronomyViewerGUI.wl` | Interactive notebook UI via `DynamicModule` (catalog dropdown + manual coordinates). |
| `interactive_sky_map.py` | `InteractiveSkyMap.wl` | Clickable all-sky map with telescope view (`DynamicModule`, `EventHandler`). |

## Running

### Command line (batch)

```bash
wolframscript -file Wolfram/HorseheadCompleteDemo.wl
wolframscript -file Wolfram/HorseheadSkyView.wl
```

Network-based scripts require outbound HTTPS to `skyview.gsfc.nasa.gov`, `skyserver.sdss.org`, or `archive.stsci.edu`.

### Notebook (interactive)

```wolfram
Get["Wolfram/AstronomyViewerGUI.wl"];
AstronomyViewerGUI`launchAstronomyViewer[]

Get["Wolfram/InteractiveSkyMap.wl"];
InteractiveSkyMap`launchInteractiveSkyMap[]
```

## Python-to-Wolfram cheat sheet

| Python | Wolfram Language |
|---|---|
| `requests.get(url, params=..)` | `URLExecute[url <> "?" <> URLQueryEncode[params], "ByteArray"]` |
| `astropy.io.fits.open(...)` | `Import[file, {"FITS", "Data"}]` (+ `"Metadata"`) |
| `numpy.linspace(a, b, n)` | `Subdivide[a, b, n-1]` |
| `numpy.meshgrid(x, y)` | `ConstantArray[x, Length[y]]` / `Transpose[ConstantArray[y, Length[x]]]` |
| `numpy.percentile(a, q)` | `Quantile[a, q/100]` |
| `numpy.random.poisson(\[Lambda])` | `RandomVariate[PoissonDistribution[\[Lambda]]]` |
| `numpy.random.normal(\[Mu], \[Sigma])` | `RandomVariate[NormalDistribution[\[Mu], \[Sigma]]]` |
| `matplotlib.imshow(data)` | `ArrayPlot[Reverse[data], ColorFunction -> ...]` |
| `matplotlib.colorbar` | `ArrayPlot[..., PlotLegends -> Automatic]` |
| `tkinter` UI | `DynamicModule[{...}, Grid[{{controls, display}}]]` |
| `mpl_connect('button_press_event', cb)` | `EventHandler[g, {"MouseClicked" :> (...)}]` |
| `json.dump(obj, f)` | `Export[file, assoc, "JSON"]` |

## Notes / Caveats

- Wolfram's built-in `Import[..., "FITS"]` handles the primary HDU image data. For full FITS/WCS support you may need the `Astronomy``` paclet or to parse header cards manually.
- Poisson photon noise is approximated with a Gaussian of matching mean and variance for speed in the large simulated image; swap in `RandomVariate[PoissonDistribution[\[Mu]]]` if you want exact statistics.
- The Tk GUIs become `DynamicModule` front-end interfaces. In batch mode (`wolframscript`) they print a reminder but do not render.

See the top-level `README.md` and `CLAUDE.md` for project context and licensing information.
