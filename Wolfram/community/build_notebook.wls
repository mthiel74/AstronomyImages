#!/usr/bin/env wolframscript
(* Build Wolfram/community/AstronomyImagesPost.nb \[LongDash] a Wolfram
   Community post that walks through the Python-to-Wolfram translation
   of the AstronomyImages toolkit, ships a dozen SkyView cutouts of
   famous targets, a 2MASS J/H/K false-colour composite, a
   multi-wavelength Horsehead panel, and a live downloader you can
   aim at any object. *)

hereDir = DirectoryName[If[$ScriptCommandLine =!= {}, $InputFileName,
                           NotebookFileName[]]];
imgDir  = FileNameJoin[{hereDir, "images"}];

imgWidth = 680;

loadImg[name_] := With[
  {p = FileNameJoin[{imgDir, name}]},
  If[FileExistsQ[p],
    Image[Import[p], ImageSize -> imgWidth],
    Style["[missing image: " <> name <> "]", Red]]];

imgCell[name_, caption_String: ""] := Cell[
  BoxData[ToBoxes[
    If[caption === "",
      loadImg[name],
      Column[{loadImg[name],
              Pane[Style[caption, Italic, FontSize -> 10,
                         FontFamily -> "Helvetica"],
                   imgWidth, Alignment -> Center]},
             Alignment -> Center, Spacings -> 0.6]]]],
  "Output", ShowCellLabel -> False,
  CellMargins -> {{0, 0}, {10, 14}}];

(* Two-column cheat-sheet grid *)
cheatGrid = Grid[{
  {Style["Python source", Bold],                 Style["Wolfram port", Bold]},
  {"api_examples.json",                          "APIExamples.wl"},
  {"horsehead_skyview.py",                       "HorseheadSkyView.wl"},
  {"horsehead_skyview_fixed.py",                 "HorseheadSkyViewFixed.wl"},
  {"horsehead_alternative.py",                   "HorseheadAlternative.wl"},
  {"horsehead_complete_demo.py",                 "HorseheadCompleteDemo.wl"},
  {"astronomy_viewer_gui.py",                    "AstronomyViewerGUI.wl"},
  {"interactive_sky_map.py",                     "InteractiveSkyMap.wl"}},
  Frame -> All, Background -> {None, {{LightGray}, White}},
  Alignment -> {Left, Center}, Spacings -> {1.2, 0.6}];

(* Licensing table *)
licenceGrid = Grid[{
  {Style["Survey", Bold], Style["Status", Bold],
   Style["Commercial use", Bold], Style["Required attribution", Bold]},
  {"2MASS J/H/K",      "Public domain (NASA/NSF)",
   "Yes",
   "Joint project of UMass and IPAC/Caltech, NASA/NSF-funded"},
  {"WISE W1\[Dash]W4", "Public domain (NASA)",
   "Yes",
   "UCLA + JPL/Caltech, NASA-funded"},
  {"GALEX",            "Public domain (NASA)",
   "Yes",
   "NASA Galaxy Evolution Explorer"},
  {"DSS / DSS2",       "Free for research & education, attribution required (STScI)",
   "Ask STScI",
   "Produced at STScI under grant NAG W-2166 from Palomar/UK Schmidt photographic data"},
  {"SDSS imagery",     "Free with acknowledgement",
   "With attribution",
   "See sdss.org/collaboration/citing-sdss"}},
  Frame -> All,
  Background -> {None, {{LightYellow},
                        {GrayLevel[0.95]}, {White}}},
  Alignment -> {Left, Top},
  Spacings -> {1.2, 0.8}];

cells = {

 Cell["From Python to Wolfram: Exploring NASA SkyView", "Title"],

 Cell["A walkthrough of the AstronomyImages toolkit port, with a \
Messier gallery and a live survey downloader", "Subtitle"],

 Cell[TextData[{
   "Marco Thiel, April 2026 \[LongDash] ",
   ButtonBox["github.com/mthiel74/AstronomyImages",
     BaseStyle -> "Hyperlink",
     ButtonData -> {URL["https://github.com/mthiel74/AstronomyImages"],
                    None}]}],
      "Text", TextAlignment -> Center, FontSlant -> Italic],

 (* ---- Abstract ---- *)

 Cell["Abstract", "Section"],

 Cell[TextData[{
   "Amateur and professional astronomy data are wildly more \
accessible than they were 20 years ago. Ten minutes of browsing the \
NASA SkyView Virtual Observatory will put half a century of \
photographic sky surveys, three near-infrared bands from 2MASS, four \
mid-infrared bands from WISE, and a pile of X-ray and radio surveys \
on your desktop, all for free, all in FITS, and \[LongDash] for the \
surveys NASA operates directly \[LongDash] all in the public domain. \
The friction is now entirely in the tooling: how do you wire the HTTP \
query, decode the FITS, stretch the pixels sensibly, overlay analysis \
widgets, and avoid reinventing a matplotlib-plus-Tkinter wheel every \
time you want to look at a new target?"}], "Text"],

 Cell[TextData[{
   "The open-source project ",
   StyleBox["AstronomyImages", "Input"],
   " is a small Python answer to that question. This notebook is its \
Wolfram-Language twin. We walk through the Python\[Dash]to\[Dash]\
Wolfram translation (",
   StyleBox["URLRead", "Input"], ", ",
   StyleBox["Import[_, \"FITS\"]", "Input"], ", ",
   StyleBox["ArrayPlot", "Input"], ", ",
   StyleBox["Manipulate", "Input"],
   "), pull a fresh 0.5\[Degree] DSS2 Red cutout of the Horsehead \
Nebula from NASA SkyView, build a 2MASS J/H/K false-colour composite \
with an asinh stretch, tour six Messier targets (M1, M16, M31, M42, \
M45, M51) across multiple colour maps, and expose a self-contained \
interactive downloader that lets you aim the whole pipeline at any \
sky coordinate. Everything in the notebook is evaluable offline \
thanks to bundled cutouts; the live cells touch the network only when \
you ask them to."}], "Text"],

 imgCell["horsehead_DSS2_Red.png",
  "Figure 1. A live 600\[Times]600 DSS2 Red cutout of the Horsehead \
Nebula (Barnard 33), centred at RA 85.25\[Degree], \
Dec \[Minus]2.46\[Degree], 0.5\[Degree] on a side. The bright star at \
the upper left is \[Sigma] Orionis; the horse-head silhouette is a \
foreground molecular cloud eating into the emission nebula IC 434 \
behind it. Rendered by Wolfram in a single ArrayPlot call. \
Image credit: Digitized Sky Survey 2 (Caltech Palomar plates), \
produced at STScI under U.S. Government grant NAG W-2166 and \
distributed through NASA SkyView; reproduced here under STScI's \
research/education terms \[LongDash] see \[Section]10."],

 (* ---- Section 1: What the toolkit does ---- *)

 Cell["1. What the original Python toolkit does", "Section"],

 Cell[TextData[{"The Python project in the parent directory contains \
five kinds of scripts:\n",
   "\[Bullet] ",
   ButtonBox["NASA SkyView",
     BaseStyle -> "Hyperlink",
     ButtonData -> {URL["https://skyview.gsfc.nasa.gov/"], None}],
   " FITS cutout downloaders (30+ surveys \[LongDash] DSS2, 2MASS, \
WISE, GALEX, RASS, NVSS\[Ellipsis]);\n",
   "\[Bullet] fallback JPEG queries against SDSS SkyServer and STScI \
DSS2 when SkyView is down;\n",
   "\[Bullet] a physics-based simulation of the Horsehead dark nebula \
\[LongDash] useful for pedagogy, and for demos in environments with \
no network;\n",
   "\[Bullet] two Tkinter GUIs (a catalog-driven viewer and a \
clickable all-sky map);\n",
   "\[Bullet] analysis scaffolding \[LongDash] histograms, \
percentile clipping, optical-depth estimates, edge detection.\n\n",
   "The Wolfram port mirrors every script as a self-contained ",
   StyleBox[".wl", "Input"],
   " package in the ", StyleBox["Wolfram/", "Input"],
   " folder. The table below shows the 1-to-1 mapping."}], "Text"],

 Cell[BoxData[ToBoxes[cheatGrid]],
  "Output", ShowCellLabel -> False,
  CellMargins -> {{0, 0}, {8, 14}}],

 (* ---- Section 2: cheat sheet ---- *)

 Cell["2. A Python\[Dash]to\[Dash]Wolfram cheat-sheet", "Section"],

 Cell[TextData[{
   "Most of the translation is a pleasant 1-to-1 substitution of \
Python idioms for Wolfram ones. A few non-obvious ones are worth \
calling out:\n",
   "\[Bullet] ",
   StyleBox["requests.get(url, params=..).content", "Input"],
   " \[Rule] ",
   StyleBox["URLRead[HTTPRequest[url <> \"?\" <> URLQueryEncode[params]]][\"BodyByteArray\"]",
     "Input"],
   "\n",
   "\[Bullet] ",
   StyleBox["astropy.io.fits.open(..)", "Input"],
   " \[Rule] ",
   StyleBox["Import[file, {\"FITS\", \"Data\"}]", "Input"],
   " (+ ",
   StyleBox["\"Metadata\"", "Input"],
   " for the header)\n",
   "\[Bullet] ",
   StyleBox["numpy.linspace(a, b, n)", "Input"],
   " \[Rule] ",
   StyleBox["Subdivide[a, b, n - 1]", "Input"],
   "\n",
   "\[Bullet] ",
   StyleBox["numpy.percentile(a, q)", "Input"],
   " \[Rule] ",
   StyleBox["Quantile[a, q/100]", "Input"],
   "\n",
   "\[Bullet] ",
   StyleBox["numpy.random.normal / poisson", "Input"],
   " \[Rule] ",
   StyleBox["RandomVariate[NormalDistribution[..]]", "Input"],
   " / ",
   StyleBox["RandomVariate[PoissonDistribution[..]]", "Input"],
   "\n",
   "\[Bullet] ",
   StyleBox["matplotlib.imshow + colorbar", "Input"],
   " \[Rule] ",
   StyleBox["ArrayPlot[Reverse[data], PlotLegends -> Automatic]",
     "Input"],
   "\n",
   "\[Bullet] ",
   StyleBox["Tkinter widget tree", "Input"],
   " \[Rule] ",
   StyleBox["DynamicModule[..] + Manipulate[..]", "Input"],
   "\n",
   "\[Bullet] ",
   StyleBox["scipy.ndimage.gaussian_filter / sobel", "Input"],
   " \[Rule] ",
   StyleBox["GaussianFilter / GradientFilter", "Input"],
   " (act directly on 2-D arrays)\n",
   "\[Bullet] ",
   StyleBox["json.dump", "Input"],
   " \[Rule] ",
   StyleBox["Export[file, assoc, \"JSON\"]", "Input"]}], "Text"],

 Cell[TextData[{
   "A subtlety worth flagging: ",
   StyleBox["URLExecute[url, \"ByteArray\"]", "Input"],
   " is a very tempting translation of ",
   StyleBox["requests.get(...).content", "Input"],
   " \[LongDash] but the second argument of ",
   StyleBox["URLExecute", "Input"],
   " is the ",
   StyleBox["Import", "Input"],
   " element, not the response encoding. When the server returns ",
   StyleBox["Content-Type: image/fits", "Input"],
   " the call therefore tries to parse the body as a FITS file and \
extract a field called \"ByteArray\", which does not exist. The \
working idiom, ",
   StyleBox["URLRead[HTTPRequest[url]][\"BodyByteArray\"]", "Input"],
   ", returns the raw bytes and plays nicely with ",
   StyleBox["Import", "Input"],
   " downstream."}], "Text"],

 (* ---- Section 3: downloading ---- *)

 Cell["3. Downloading a real image from NASA SkyView", "Section"],

 Cell[TextData[{"Every parameter SkyView needs is packed into the \
query string: where on the sky (",
   StyleBox["Position", "Input"],
   ", in equatorial degrees), which survey, how many pixels on a \
side, how big a field (",
   StyleBox["Size", "Input"],
   ", in degrees), and what to return. The cell below builds the URL \
and fetches the bytes."}], "Text"],

 Cell[BoxData["\
skyViewURL[ra_, dec_, survey_, fov_, pixels_] :=\n\
  \"https://skyview.gsfc.nasa.gov/cgi-bin/images?\" <>\n\
    URLQueryEncode[<|\n\
      \"Position\" -> ToString[ra] <> \",\" <> ToString[dec],\n\
      \"Survey\"   -> survey,\n\
      \"Pixels\"   -> ToString[pixels] <> \",\" <> ToString[pixels],\n\
      \"Size\"     -> ToString[fov],\n\
      \"Return\"   -> \"FITS\"\n\
    |>];\n\
\n\
skyViewBytes[ra_, dec_, survey_, fov_: 0.5, pixels_: 600] :=\n\
  Module[{resp = URLRead[HTTPRequest[skyViewURL[ra, dec, survey, fov, pixels]]]},\n\
    If[resp === $Failed || resp[\"StatusCode\"] =!= 200,\n\
      $Failed, resp[\"BodyByteArray\"]]];"], "Input"],

 Cell[TextData[{"The notebook ships with eighteen pre-downloaded \
cutouts in ",
   StyleBox["images/", "Input"],
   " so every cell below can be re-run offline. The next cell loads \
the DSS2 Red cutout of the Horsehead and turns it into a numeric \
matrix. ",
   StyleBox["Import[_, {\"FITS\", \"Data\"}]", "Input"],
   " returns a one-entry Association (",
   StyleBox["<|1 -> array|>", "Input"],
   ") for primary-HDU files and a longer Association for multi-HDU \
files; ",
   StyleBox["First[Values[..]]", "Input"],
   " grabs the primary HDU either way."}], "Text"],

 Cell[BoxData["\
loadFITSArray[file_] :=\n\
  With[{raw = Import[file, {\"FITS\", \"Data\"}]},\n\
    Which[\n\
      AssociationQ[raw],                First[Values[raw]],\n\
      ListQ[raw] && Depth[raw] >= 4,    First[raw],\n\
      True,                             raw]];\n\
\n\
horseheadFile = FileNameJoin[{NotebookDirectory[],\n\
                              \"images\", \"horsehead_DSS2_Red.fits\"}];\n\
data = loadFITSArray[horseheadFile];\n\
Dimensions[data]"], "Input"],

 Cell[TextData[{"A quick percentile-clipped rendering with ",
   StyleBox["ArrayPlot", "Input"],
   " gives a publication-quality image in one line. ",
   StyleBox["Reverse", "Input"],
   " flips the array so the FITS origin (lower-left) ends up where \
astronomers expect it."}], "Text"],

 Cell[BoxData["\
quickView[d_, cf_: \"GrayTones\"] :=\n\
  Module[{flat = Flatten[d], lo, hi},\n\
    {lo, hi} = Quantile[flat, {0.01, 0.995}];\n\
    ArrayPlot[Reverse[d],\n\
      ColorFunction -> ColorData[cf],\n\
      PlotRange -> {lo, hi},\n\
      ImageSize -> 480, Frame -> False]];\n\
\n\
quickView[data]"], "Input"],

 (* ---- Section 4: multi-wavelength ---- *)

 Cell["4. One object, four wavelengths", "Section"],

 Cell[TextData[{"SkyView stitches together data from many surveys, so \
with one coordinate and five lines of Wolfram you can look at the \
same patch of sky in completely different wavelengths. The Horsehead \
makes the point vividly: it is a ",
   StyleBox["dark", FontSlant -> Italic],
   " nebula \[LongDash] an opaque column of dust that absorbs the \
red H\[Alpha] light of IC 434 behind it \[LongDash] so it jumps out \
as a silhouette at visible wavelengths. In the near-infrared the same \
dust is \[TildeTilde]10 times less absorbing; in the mid-infrared \
the dust is ",
   StyleBox["emitting", FontSlant -> Italic],
   ", so the nebula flips from dark feature to bright feature as we \
slide from 656 nm to 3.4 \[Mu]m."}], "Text"],

 imgCell["horsehead_multiwavelength.png",
  "Figure 2. The Horsehead Nebula in four bands spanning a factor of \
~5000 in wavelength. DSS2 Red (656 nm) shows the classic silhouette \
absorbing H\[Alpha]. 2MASS J (1.25 \[Mu]m) is only 2\[Times] redward \
but the extinction is already half its optical value \[LongDash] \
field stars punch through. At 2MASS K (2.17 \[Mu]m) the nebula is \
almost transparent. WISE W1 (3.4 \[Mu]m, shown in SunsetColors) \
reveals warm dust emission from the nebula itself, plus the \
broadband glow of PAH molecules along the IC 434 ridge. Image \
credits: DSS2 / STScI (research-education use, see \[Section]10); \
2MASS J, K \[Dash] IPAC/Caltech, UMass (public domain, NASA/NSF); \
WISE W1 \[Dash] NASA / JPL-Caltech / UCLA (public domain)."],

 Cell[TextData[{
   "The dust extinction curve is roughly ",
   Cell[BoxData[FormBox[RowBox[{
     SubscriptBox["A", "\[Lambda]"], "\[Proportional]",
     SuperscriptBox["\[Lambda]",
       RowBox[{"-", StyleBox["\[Beta]", FontSlant -> "Italic"]}]]}],
     TraditionalForm]]],
   " with \[Beta] \[TildeTilde] 1.5 in the optical/near-IR. Plugging \
in the Horsehead's optical-depth-of-order-unity at 550 nm gives \
\[Tau] \[TildeTilde] 0.1 at 2 \[Mu]m \[LongDash] consistent with the \
2MASS-K panel above, where the dust is barely visible."}], "Text"],

 (* ---- Section 5: RGB composite ---- *)

 Cell["5. A 2MASS J/H/K false-colour composite", "Section"],

 Cell[TextData[{"The three 2MASS bands (J at 1.25 \[Mu]m, H at 1.65, \
K at 2.17) were taken with the same camera at the same time through \
a single filter wheel, so they share pixel scale, pointing, and \
seeing \[LongDash] a rare luxury in astronomy. That makes them the \
textbook ingredients for a false-colour composite. We map J to blue, \
H to green, K to red; the colour of a pixel then encodes its \
spectral shape, with hot dust and reddened objects leaning K and \
therefore looking red in the final image."}], "Text"],

 Cell[TextData[{"Raw telescope counts span many orders of magnitude \
\[LongDash] a star can be 10,000\[Times] brighter than the sky \
background it sits on. A linear stretch wastes most of the colour \
space on a handful of stars. The solution is an ",
   StyleBox["asinh", "Input"],
   " (inverse hyperbolic sine) stretch, which is linear near zero \
and logarithmic at large values: it preserves the shape of faint \
features while still fitting bright stars into the [0,1] interval."}],
      "Text"],

 Cell[BoxData["\
asinhStretch[data_, qLo_: 0.5, qHi_: 0.998, softening_: 0.02] :=\n\
  Module[{flat = Flatten[data], lo, hi, norm, s, sMin, sMax},\n\
    {lo, hi} = Quantile[flat, {qLo, qHi}];\n\
    norm = Clip[(data - lo)/(hi - lo + 10.^-12), {0, 1}];\n\
    s = ArcSinh[norm/softening];\n\
    {sMin, sMax} = MinMax[Flatten[s]];\n\
    Clip[(s - sMin)/(sMax - sMin + 10.^-12), {0, 1}]];\n\
\n\
{jD, hD, kD} = loadFITSArray /@ {\n\
  FileNameJoin[{NotebookDirectory[], \"images\", \"horsehead_2MASS_J.fits\"}],\n\
  FileNameJoin[{NotebookDirectory[], \"images\", \"horsehead_2MASS_H.fits\"}],\n\
  FileNameJoin[{NotebookDirectory[], \"images\", \"horsehead_2MASS_K.fits\"}]};\n\
\n\
rgb = Transpose[Map[Reverse,\n\
        {asinhStretch[kD], asinhStretch[hD], asinhStretch[jD]}],\n\
       {3, 1, 2}];\n\
Image[rgb, Interleaving -> True, ImageSize -> 520]"], "Input"],

 imgCell["horsehead_RGB_JHK.png",
  "Figure 3. Near-infrared false colour: 2MASS K (2.17 \[Mu]m) in red, \
H (1.65 \[Mu]m) in green, J (1.25 \[Mu]m) in blue. Reddened stars \
behind the dust column glow orange; the nebular emission of IC 434 \
shines blue-white. An asinh stretch keeps stars from saturating \
while preserving low-surface-brightness structure. Image credit: \
Two Micron All Sky Survey (IPAC/Caltech + UMass, NASA/NSF-funded), \
public domain."],

 (* ---- Section 6: rudimentary analysis ---- *)

 Cell["6. Rudimentary image analysis", "Section"],

 Cell[TextData[{"Every tool you need for a first-pass numerical look \
at a FITS image is in the top-level ",
   StyleBox["System`", "Input"],
   " context. No toolbox imports."}], "Text"],

 Cell["6.1 Summary statistics", "Subsection"],

 Cell[BoxData["flat = Flatten[data];"], "Input"],

 Cell[BoxData["\
<|\"Min\"    -> Min[flat],\n\
  \"Max\"    -> Max[flat],\n\
  \"Mean\"   -> N @ Mean[flat],\n\
  \"Median\" -> N @ Median[flat],\n\
  \"StdDev\" -> N @ StandardDeviation[flat]|>"], "Input"],

 Cell["6.2 Intensity histogram", "Subsection"],

 Cell[BoxData["\
Histogram[flat, 80,\n\
  ChartStyle -> RGBColor[0.27, 0.51, 0.71],\n\
  Frame -> True,\n\
  FrameLabel -> {\"Pixel intensity (ADU)\", \"Count\"},\n\
  PlotLabel -> \"Intensity distribution\",\n\
  ImageSize -> 520]"], "Input"],

 Cell["6.3 Dark-nebula segmentation and optical depth", "Subsection"],

 Cell[TextData[{"The Horsehead's dust column attenuates background \
starlight by a factor ",
   Cell[BoxData[FormBox[RowBox[{
     SuperscriptBox["e", RowBox[{"-", "\[Tau]"}]]}],
     TraditionalForm]]],
   " where \[Tau] is the optical depth. A crude segmentation flags \
every pixel below the 30th percentile; comparing the mean intensity \
inside and outside the mask gives a direct estimate of \[Tau] via \
the Beer-Lambert law, ",
   Cell[BoxData[FormBox[RowBox[{SubscriptBox["I", "dark"], "=",
     SubscriptBox["I", "bg"],
     SuperscriptBox["e", RowBox[{"-", "\[Tau]"}]]}],
     TraditionalForm]]],
   ". For DSS2 Red we expect \[Tau] \[TildeTilde] 0.5 for the dense \
part of the column."}], "Text"],

 Cell[BoxData["threshold = Quantile[flat, 0.30]"], "Input"],

 Cell[BoxData["darkMask = UnitStep[threshold - data];"], "Input"],

 Cell[BoxData["\
{quickView[data],\n\
 ArrayPlot[Reverse[darkMask],\n\
   ColorRules -> {0 -> White, 1 -> Red},\n\
   ImageSize -> 480, Frame -> False]}"], "Input"],

 Cell[BoxData["\
bgMean  = Mean[Pick[flat, Flatten[darkMask], 0]];\n\
nebMean = Mean[Pick[flat, Flatten[darkMask], 1]];\n\
tauEst  = -Log[nebMean/bgMean];\n\
<|\"bg (ADU)\"        -> N @ bgMean,\n\
  \"dark (ADU)\"      -> N @ nebMean,\n\
  \"optical depth\"   -> N @ tauEst,\n\
  \"extinction (mag)\" -> N[1.086 tauEst]|>"], "Input"],

 imgCell["horsehead_analysis_panel.png",
  "Figure 4. Four-panel summary of the DSS2 Red cutout: raw \
percentile-clipped image (upper left), log10-scaled version revealing \
faint IC 434 emission (upper right), dark-nebula segmentation at the \
30th percentile (lower left), and intensity histogram with the \
threshold marked (lower right). Image credit: DSS2 / STScI \
(research-education use, see \[Section]10)."],

 Cell["6.4 Smoothing, log scaling, edge detection", "Subsection"],

 Cell[TextData[{"Wolfram's image-processing operators all take plain \
2-D arrays \[LongDash] no conversion needed. ",
   StyleBox["GaussianFilter", "Input"],
   " smooths at a given pixel scale; ",
   StyleBox["GradientFilter", "Input"],
   " returns |\[Del] ",
   StyleBox["I", FontSlant -> Italic],
   "|, which lights up edges."}], "Text"],

 Cell[BoxData["\
Grid[{{\n\
  Labeled[ArrayPlot[Reverse[Log10[data + 1]],\n\
      ColorFunction -> \"SunsetColors\",\n\
      ImageSize -> 320, Frame -> False], \"log10 scaled\"],\n\
  Labeled[ArrayPlot[Reverse[GaussianFilter[data, 3]],\n\
      ColorFunction -> \"GrayTones\",\n\
      ImageSize -> 320, Frame -> False], \"Gaussian \[Sigma]=3\"],\n\
  Labeled[ArrayPlot[Reverse[GradientFilter[data, 2]],\n\
      ColorFunction -> \"TemperatureMap\",\n\
      ImageSize -> 320, Frame -> False], \"Gradient (edge)\"]}}]"],
      "Input"],

 (* ---- Section 7: Messier gallery ---- *)

 Cell["7. A Messier gallery", "Section"],

 Cell[TextData[{"To convince you the pipeline isn't Horsehead-only, \
below is a six-target tour through the Messier catalogue, all \
fetched from NASA SkyView in DSS2 Red at the same resolution, each \
rendered with a different colour map to show how much the ",
   StyleBox["same", FontSlant -> Italic],
   " raw data can change character with a one-line style swap. \
Every object has its own challenge \[LongDash] M31 is so large it \
needs 3\[Degree] on a side, M1 is compact at ~5\[CapitalAHat]', \
Pleiades is dominated by a handful of bright B-type stars that \
saturate the central pixels on the original plates."}], "Text"],

 imgCell["messier_gallery.png",
  "Figure 5. Six famous Messier targets fetched from NASA SkyView, \
all in DSS2 Red at the same resolution, each rendered with a \
different colour map. From upper-left, reading rows: M42 \
\[LongDash] the Orion Nebula, ~1344 ly away, the nearest active \
stellar nursery. M45 \[LongDash] the Pleiades, a ~110 Myr-old open \
cluster of B stars in Taurus. M31 \[LongDash] the Andromeda Galaxy, \
our nearest spiral at 2.5 Mly, on a collision course with the Milky \
Way. M16 \[LongDash] the Eagle Nebula, home of Hubble's Pillars of \
Creation. M51 \[LongDash] the Whirlpool Galaxy interacting with NGC \
5195 at 23 Mly. M1 \[LongDash] the Crab Nebula, the expanding \
remnant of a supernova seen by Chinese astronomers in 1054 AD. \
Image credit: Digitized Sky Survey 2 (Caltech Palomar + UK Schmidt \
plates), produced at STScI under grant NAG W-2166 and distributed \
through NASA SkyView; reproduced under STScI's research-education \
terms (see \[Section]10)."],

 Cell[TextData[{"For individual high-resolution plates of each \
target, see ",
   StyleBox["images/gallery_*.png", "Input"],
   " in the notebook directory, or re-run ",
   StyleBox["make_assets.wls", "Input"],
   " to regenerate them from scratch."}], "Text"],

 (* ---- Section 8: interactive explorer ---- *)

 Cell["8. Interactive explorers", "Section"],

 Cell["8.1 Offline Manipulate over the bundled surveys", "Subsection"],

 Cell[TextData[{"The ",
   StyleBox["Manipulate", "Input"],
   " below turns Section 6's analysis into a live tool: scrub \
through three bundled cutouts, change the colour map, and slide the \
dark-threshold percentile up or down while the segmentation and \
histogram update in real time. No network traffic \[LongDash] \
everything reads from ",
   StyleBox["images/", "Input"], "."}], "Text"],

 Cell[BoxData["\
surveyFiles = <|\n\
  \"DSS2 Red\"   -> FileNameJoin[{NotebookDirectory[], \"images\",\n\
                                    \"horsehead_DSS2_Red.fits\"}],\n\
  \"2MASS-J\"    -> FileNameJoin[{NotebookDirectory[], \"images\",\n\
                                    \"horsehead_2MASS_J.fits\"}],\n\
  \"Simulated\"  -> FileNameJoin[{NotebookDirectory[], \"images\",\n\
                                    \"horsehead_simulated.mx\"}]|>;\n\
\n\
loadImage[f_] :=\n\
  If[StringEndsQ[f, \".mx\"], Import[f], loadFITSArray[f]];"], "Input"],

 Cell[BoxData["\
Manipulate[\n\
  Module[{img = loadImage[surveyFiles[sel]], flatImg, lo, hi, thr, mask},\n\
    flatImg = Flatten[img];\n\
    {lo, hi} = Quantile[flatImg, {0.01, 0.995}];\n\
    thr  = Quantile[flatImg, pct/100.];\n\
    mask = UnitStep[thr - img];\n\
    Grid[{{\n\
      ArrayPlot[Reverse[img],\n\
        ColorFunction -> ColorData[cf],\n\
        PlotRange -> {lo, hi},\n\
        ImageSize -> 380, Frame -> False],\n\
      Show[\n\
        ArrayPlot[Reverse[img],\n\
          ColorFunction -> ColorData[cf],\n\
          PlotRange -> {lo, hi},\n\
          ImageSize -> 380, Frame -> False],\n\
        ArrayPlot[Reverse[mask],\n\
          ColorRules -> {0 -> Opacity[0],\n\
                         1 -> RGBColor[1, 0, 0, 0.45]},\n\
          ImageSize -> 380, Frame -> False]],\n\
      Histogram[flatImg, 60,\n\
        ChartStyle -> RGBColor[0.27, 0.51, 0.71],\n\
        GridLines -> {{{thr, Directive[Red, Dashed]}}, None},\n\
        Frame -> True, ImageSize -> 380]}}]],\n\
  Style[\"Survey\", Bold],\n\
  {{sel, \"DSS2 Red\"}, {\"DSS2 Red\", \"2MASS-J\", \"Simulated\"}},\n\
  Style[\"Colour map\", Bold],\n\
  {{cf, \"GrayTones\"},\n\
    {\"GrayTones\", \"SunsetColors\", \"TemperatureMap\",\n\
     \"Rainbow\", \"DeepSeaColors\"}},\n\
  Style[\"Dark-threshold percentile\", Bold],\n\
  {{pct, 30}, 5, 70, 1}]"], "Input"],

 Cell["8.2 Live downloader \[LongDash] any object, any survey",
      "Subsection"],

 Cell[TextData[{"The next cell is a small interactive viewer that \
actually hits the network when you press the button. Pick a target \
from the catalog (or type in your own RA/Dec), choose a survey and a \
field of view, and the cell will fetch the FITS from NASA SkyView, \
decode it, and render it with your chosen colour map. Everything \
lives inside one ",
   StyleBox["DynamicModule", "Input"],
   " \[LongDash] no global state leaks out. The catalog is just an ",
   StyleBox["Association", "Input"],
   " so you can extend it to any target you like. The figure below \
shows what the evaluated cell looks like after a fresh download of \
M42 in DSS2 Red."}], "Text"],

 imgCell["gui_preview.png",
  "Figure 7. The live downloader after pressing Download & plot with \
M42 Orion Nebula / DSS2 Red / 1\[Degree] FOV / 600 px / SunsetColors. \
Every control is a built-in Wolfram widget \[LongDash] a catalog \
picker, RA/Dec input fields, an FOV slider, a survey picker, a \
pixel-size picker, and a colour-map picker \[LongDash] wired up \
inside a single DynamicModule. The status line under the buttons \
reports success, image size, and bytes received. Image credit: \
DSS2 / STScI (research-education use, see \[Section]10)."],

 Cell[BoxData["\
catalog = <|\n\
  \"Horsehead Nebula\"  -> <|\"ra\" -> 85.2479,  \"dec\" -> -2.4583,  \"fov\" -> 0.5|>,\n\
  \"M42 Orion Nebula\"  -> <|\"ra\" -> 83.8213,  \"dec\" -> -5.3903,  \"fov\" -> 1.0|>,\n\
  \"M45 Pleiades\"      -> <|\"ra\" -> 56.870,   \"dec\" -> 24.105,   \"fov\" -> 2.5|>,\n\
  \"M31 Andromeda\"     -> <|\"ra\" -> 10.685,   \"dec\" -> 41.269,   \"fov\" -> 3.0|>,\n\
  \"M16 Eagle Nebula\"  -> <|\"ra\" -> 274.700,  \"dec\" -> -13.807,  \"fov\" -> 0.5|>,\n\
  \"M51 Whirlpool\"     -> <|\"ra\" -> 202.469,  \"dec\" -> 47.196,   \"fov\" -> 0.2|>,\n\
  \"M1 Crab Nebula\"    -> <|\"ra\" -> 83.633,   \"dec\" -> 22.014,   \"fov\" -> 0.15|>,\n\
  \"M104 Sombrero\"     -> <|\"ra\" -> 189.998,  \"dec\" -> -11.623,  \"fov\" -> 0.15|>,\n\
  \"NGC 2024 Flame\"    -> <|\"ra\" -> 85.425,   \"dec\" -> -1.900,   \"fov\" -> 0.4|>,\n\
  \"Veil Nebula\"       -> <|\"ra\" -> 312.75,   \"dec\" -> 30.72,    \"fov\" -> 3.0|>,\n\
  \"(custom)\"          -> <|\"ra\" -> 0.,       \"dec\" -> 0.,       \"fov\" -> 0.5|>\n\
  |>;\n\
\n\
surveys = {\"DSS2 Red\", \"DSS2 Blue\", \"DSS2 IR\",\n\
           \"2MASS-J\", \"2MASS-H\", \"2MASS-K\",\n\
           \"WISE 3.4\", \"WISE 4.6\", \"WISE 12\", \"WISE 22\",\n\
           \"GALEX Near UV\"};"], "Input"],

 Cell[BoxData["\
DynamicModule[{\n\
    pick = \"Horsehead Nebula\",\n\
    ra = 85.2479, dec = -2.4583, fov = 0.5,\n\
    survey = \"DSS2 Red\", pixels = 600, cmap = \"GrayTones\",\n\
    status = \"Ready.\",\n\
    data = None, lastURL = \"\"},\n\
  (* Keep RA/Dec/FOV in sync with the catalog picker *)\n\
  Row[{\n\
    Panel[Column[{\n\
      Style[\"Live SkyView downloader\", Bold, 14],\n\
      Row[{\"Target: \", PopupMenu[Dynamic[pick,\n\
         Function[v, pick = v;\n\
           If[v =!= \"(custom)\",\n\
             ra  = catalog[v][\"ra\"];\n\
             dec = catalog[v][\"dec\"];\n\
             fov = catalog[v][\"fov\"]]]],\n\
        Keys[catalog], FieldSize -> 22]}],\n\
      Row[{\"RA (\[Degree]): \",  InputField[Dynamic[ra],  Number, FieldSize -> 10]}],\n\
      Row[{\"Dec (\[Degree]): \", InputField[Dynamic[dec], Number, FieldSize -> 10]}],\n\
      Row[{\"FOV (\[Degree]): \", Slider[Dynamic[fov], {0.02, 5.0, 0.02}],\n\
         \" \", Dynamic[NumberForm[fov, {3, 2}]]}],\n\
      Row[{\"Survey: \", PopupMenu[Dynamic[survey], surveys, FieldSize -> 22]}],\n\
      Row[{\"Pixels: \", PopupMenu[Dynamic[pixels], {300, 500, 600, 800, 1000}]}],\n\
      Row[{\"Colour map: \", PopupMenu[Dynamic[cmap],\n\
        {\"GrayTones\", \"SunsetColors\", \"TemperatureMap\", \"Rainbow\",\n\
         \"DeepSeaColors\", \"AuroraColors\", \"SolarColors\"}]}],\n\
      Row[{\n\
        Button[Style[\"Download & plot\", Bold, White],\n\
          status = \"Contacting SkyView\[Ellipsis]\";\n\
          Module[{bytes, tmp, arr},\n\
            bytes = skyViewBytes[ra, dec, survey, fov, pixels];\n\
            lastURL = skyViewURL[ra, dec, survey, fov, pixels];\n\
            If[!ByteArrayQ[bytes] || Length[bytes] < 8000,\n\
              status = \"\[Times] SkyView returned no data (survey \
may not cover this coordinate).\"; data = None,\n\
              tmp = CreateFile[];\n\
              BinaryWrite[tmp, bytes]; Close[tmp];\n\
              arr = loadFITSArray[tmp];\n\
              DeleteFile[tmp];\n\
              If[MatrixQ[arr, NumericQ],\n\
                data = arr;\n\
                status = \"\[Checkmark] \" <> ToString[Dimensions[arr]] <>\n\
                         \" image from \" <> survey <> \" (\" <>\n\
                         ToString[Length[bytes]] <> \" bytes)\",\n\
                data = None;\n\
                status = \"\[Times] Could not decode FITS.\"]]],\n\
          Background -> RGBColor[0.12, 0.45, 0.80],\n\
          ImageSize -> 200],\n\
        Spacer[6],\n\
        Button[\"Save FITS\",\n\
          If[data =!= None,\n\
            Module[{p = SystemDialogInput[\"FileSave\", \"cutout.fits\"]},\n\
              If[StringQ[p],\n\
                BinaryWrite[p, skyViewBytes[ra, dec, survey, fov, pixels]];\n\
                Close[p];\n\
                status = \"Saved to \" <> p]]],\n\
          ImageSize -> 120]\n\
      }],\n\
      Style[Dynamic[status], Italic, FontSize -> 11]\n\
    }, Spacings -> 0.8], ImageSize -> 400],\n\
    Spacer[10],\n\
    Dynamic[\n\
      If[data === None,\n\
        Framed[Style[\n\
          \"Choose a target or enter coordinates, then click \
\\\"Download & plot\\\".\",\n\
          14, Gray], ImageSize -> {540, 540}, Alignment -> Center],\n\
        Module[{lo, hi},\n\
          {lo, hi} = Quantile[Flatten[data], {0.01, 0.995}];\n\
          ArrayPlot[Reverse[data],\n\
            ColorFunction -> ColorData[cmap],\n\
            PlotRange -> {lo, hi},\n\
            Frame -> True,\n\
            PlotLabel -> Style[pick <> \" \[LongDash] \" <> survey,\n\
                               Bold, 12],\n\
            FrameLabel -> {\"X (pixels)\", \"Y (pixels)\"},\n\
            ImageSize -> 540]]]]\n\
  }]]"], "Input"],

 Cell[TextData[{"The full all-sky clickable viewer lives in ",
   StyleBox["Wolfram/InteractiveSkyMap.wl", "Input"],
   ". Load it in a notebook and evaluate ",
   StyleBox["InteractiveSkyMap`launchInteractiveSkyMap[]", "Input"],
   " for a clickable Mercator with constellation overlays and the \
same downloader wired in."}], "Text"],

 (* ---- Section 9: simulated ---- *)

 Cell["9. When the network is down: a physics-based simulation",
      "Section"],

 Cell[TextData[{"Teaching environments often ban outbound network \
access, so the toolkit includes a physics-based simulation of the \
Horsehead that uses only ",
   StyleBox["RandomVariate", "Input"],
   ", ",
   StyleBox["Sin", "Input"],
   ", and ",
   StyleBox["Exp", "Input"],
   ". A radial Gaussian mimics the diffuse IC 434 emission; five \
overlaid sinusoids add large-scale structure at random angles and \
wavenumbers; two rotated anisotropic Gaussians play the role of the \
dust column, combined into an optical depth ",
   Cell[BoxData[FormBox[RowBox[{"\[Tau]", "(",
     StyleBox["x", FontSlant -> "Italic"], ",",
     StyleBox["y", FontSlant -> "Italic"], ")"}], TraditionalForm]]],
   "; the final image is the background times ",
   Cell[BoxData[FormBox[SuperscriptBox["e", RowBox[{"-", "\[Tau]"}]],
     TraditionalForm]]],
   ", Poisson-sampled and blurred with a small additive Gaussian to \
mimic read noise. The resulting array plugs into the same analysis \
cells as a real DSS cutout."}], "Text"],

 imgCell["horsehead_simulated.png",
  "Figure 6. The built-in simulation, identical in interface to a \
real SkyView image. Useful for teaching, demos, and any environment \
where network calls are blocked."],

 (* ---- Section 10: data rights ---- *)

 Cell["10. Data rights and attribution", "Section"],

 Cell[TextData[{"Every image in this notebook was fetched from NASA \
SkyView, which aggregates dozens of sky surveys with a range of \
licensing terms. The short version: the NASA-operated surveys (",
   StyleBox["2MASS, WISE, GALEX, ROSAT All-Sky Survey, NVSS", "Input"],
   ") are unambiguously ",
   StyleBox["public domain", FontWeight -> Bold],
   " under 17 U.S.C. \[Section] 105; photographic sky surveys \
distributed through STScI (",
   StyleBox["DSS, DSS2", "Input"],
   ") are free for research and education with a required \
attribution sentence. The table below summarises what you can do \
with each."}], "Text"],

 Cell[BoxData[ToBoxes[licenceGrid]],
  "Output", ShowCellLabel -> False,
  CellMargins -> {{0, 0}, {8, 14}}],

 Cell[TextData[{StyleBox["Required attribution for DSS imagery",
     FontWeight -> Bold],
   " (reproduced verbatim from ",
   ButtonBox["archive.stsci.edu/dss/acknowledging.html",
     BaseStyle -> "Hyperlink",
     ButtonData -> {URL["https://archive.stsci.edu/dss/acknowledging.html"],
                    None}],
   "):\n\n\[OpenCurlyDoubleQuote]The Digitized Sky Surveys were \
produced at the Space Telescope Science Institute under U.S. \
Government grant NAG W-2166. The images of these surveys are based \
on photographic data obtained using the Oschin Schmidt Telescope on \
Palomar Mountain and the UK Schmidt Telescope.\[CloseCurlyDoubleQuote]"}],
      "Text"],

 Cell[TextData[{StyleBox["Recommended attribution for 2MASS",
     FontWeight -> Bold],
   " (reproduced verbatim from ",
   ButtonBox["irsa.ipac.caltech.edu",
     BaseStyle -> "Hyperlink",
     ButtonData -> {URL["https://irsa.ipac.caltech.edu/Missions/2mass.html"],
                    None}],
   "):\n\n\[OpenCurlyDoubleQuote]This publication makes use of data \
products from the Two Micron All Sky Survey, which is a joint project \
of the University of Massachusetts and the Infrared Processing and \
Analysis Center / California Institute of Technology, funded by the \
National Aeronautics and Space Administration and the National \
Science Foundation.\[CloseCurlyDoubleQuote]"}], "Text"],

 Cell[TextData[{StyleBox["Wolfram Community compatibility. ",
     FontWeight -> Bold],
   "Posts on community.wolfram.com are redistributed by Wolfram \
under CC-BY-SA 4.0 (see the ",
   ButtonBox["Community Terms of Use",
     BaseStyle -> "Hyperlink",
     ButtonData -> {URL["https://community.wolfram.com/tou"], None}],
   "). Public-domain imagery (2MASS, WISE, GALEX, \[Ellipsis]) flows \
through that licence without friction. DSS imagery is distributed \
under STScI's standard research-and-education terms; including the \
attribution sentence above next to every DSS figure is the \
conventional way to discharge that obligation in a Community post. \
This notebook does so by (a) captioning every DSS figure with the \
STScI credit and (b) reproducing the full attribution text above."}],
      "Text"],

 (* ---- Section 11: takeaways ---- *)

 Cell["11. Takeaways", "Section"],

 Cell[TextData[{
   "\[Bullet] ",
   StyleBox["URLRead + HTTPRequest + URLQueryEncode", "Input"],
   " replace the Python ",
   StyleBox["requests", "Input"], " / ",
   StyleBox["urllib", "Input"],
   " dance in two lines.\n",
   "\[Bullet] ",
   StyleBox["Import[_, {\"FITS\", \"Data\"}]", "Input"],
   " reads SkyView's raw bytes with no ",
   StyleBox["BytesIO", "Input"],
   " shim; the only quirk is that it returns a single-entry \
Association for primary-HDU files.\n",
   "\[Bullet] ",
   StyleBox["ArrayPlot", "Input"], ", ",
   StyleBox["Histogram", "Input"], ", ",
   StyleBox["GradientFilter", "Input"],
   " and ",
   StyleBox["GaussianFilter", "Input"],
   " cover most of matplotlib + scipy.ndimage for this use case, \
all operating directly on numeric matrices.\n",
   "\[Bullet] ",
   StyleBox["DynamicModule", "Input"], " + ",
   StyleBox["PopupMenu", "Input"], " + ",
   StyleBox["Button", "Input"],
   " turns a 30-line Tkinter GUI into a 30-line notebook cell that \
works on any platform with no installation beyond Wolfram itself.\n",
   "\[Bullet] An ",
   StyleBox["asinh", "Input"],
   " stretch is the textbook recipe for turning raw telescope counts \
into a viewable RGB composite; three ",
   StyleBox["Reverse[..]", "Input"],
   " channels and a ",
   StyleBox["Transpose[_, {3, 1, 2}]", "Input"],
   " produce the final image in one line.\n",
   "\[Bullet] Images shipped with the notebook are NASA SkyView \
public-domain (2MASS, WISE) or STScI-distributed research-and-\
education data (DSS) \[LongDash] see \[Section]10."}], "Text"],

 Cell["12. References and further reading", "Section"],

 Cell[TextData[{
   "\[Bullet] NASA SkyView: ",
   ButtonBox["skyview.gsfc.nasa.gov",
     BaseStyle -> "Hyperlink",
     ButtonData -> {URL["https://skyview.gsfc.nasa.gov/"], None}],
   "\n",
   "\[Bullet] DSS acknowledgement: ",
   ButtonBox["archive.stsci.edu/dss/acknowledging.html",
     BaseStyle -> "Hyperlink",
     ButtonData -> {URL["https://archive.stsci.edu/dss/acknowledging.html"],
                    None}],
   "\n",
   "\[Bullet] 2MASS: ",
   ButtonBox["irsa.ipac.caltech.edu/Missions/2mass.html",
     BaseStyle -> "Hyperlink",
     ButtonData -> {URL["https://irsa.ipac.caltech.edu/Missions/2mass.html"],
                    None}],
   "; canonical citation Skrutskie et al. (2006), ",
   StyleBox["AJ", FontSlant -> Italic], " 131, 1163.\n",
   "\[Bullet] WISE: ",
   ButtonBox["wise2.ipac.caltech.edu",
     BaseStyle -> "Hyperlink",
     ButtonData -> {URL["https://wise2.ipac.caltech.edu/docs/release/allsky/"],
                    None}],
   "\n",
   "\[Bullet] Lupton et al. 2004 (arXiv:astro-ph/0312483) \[LongDash] \
the asinh-stretch colour-composite recipe used in \[Section]5.\n",
   "\[Bullet] This repository: ",
   ButtonBox["github.com/mthiel74/AstronomyImages",
     BaseStyle -> "Hyperlink",
     ButtonData -> {URL["https://github.com/mthiel74/AstronomyImages"],
                    None}],
   " \[LongDash] MIT-licensed code, bundled public-domain and \
research/education imagery."}], "Text"]
};

nb = Notebook[cells,
  WindowSize -> {1100, 800},
  WindowMargins -> {{Automatic, 0}, {Automatic, 0}},
  FrontEndVersion -> "15.0 for Mac OS X ARM (64-bit) (January 2026)",
  StyleDefinitions -> "Default.nb",
  Saveable -> True];

outNb  = FileNameJoin[{hereDir, "AstronomyImagesPost.nb"}];
outPdf = FileNameJoin[{hereDir, "AstronomyImagesPost.pdf"}];

Print["writing notebook -> ", outNb];
Export[outNb, nb];

Print["writing pdf      -> ", outPdf];
Export[outPdf, nb];

Print["done; nb=", FileByteCount[outNb],
      "  pdf=", FileByteCount[outPdf]];
