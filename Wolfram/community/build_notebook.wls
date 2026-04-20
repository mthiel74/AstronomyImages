#!/usr/bin/env wolframscript
(* Build Wolfram/community/AstronomyImagesPost.nb — a Wolfram Community
   post that walks through the Python-to-Wolfram translation of the
   AstronomyImages toolkit, downloads a real DSS2 Red cutout of the
   Horsehead Nebula from NASA SkyView, and exposes a Manipulate-based
   interactive explorer.

   Usage:   wolframscript -file Wolfram/community/build_notebook.wls
   Output:  Wolfram/community/AstronomyImagesPost.nb
            Wolfram/community/AstronomyImagesPost.pdf
*)

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
      Column[{
        loadImg[name],
        Pane[Style[caption, Italic, FontSize -> 10,
                   FontFamily -> "Helvetica"],
             imgWidth, Alignment -> Center]},
        Alignment -> Center, Spacings -> 0.6]]]],
  "Output", ShowCellLabel -> False,
  CellMargins -> {{0, 0}, {10, 14}}];

(* A small 2-column conversion table "Python source" -> "Wolfram port". *)
cheatGrid = Grid[{
  {Style["Python source", Bold],                 Style["Wolfram port", Bold]},
  {"api_examples.json",                          "APIExamples.wl"},
  {"horsehead_skyview.py",                       "HorseheadSkyView.wl"},
  {"horsehead_skyview_fixed.py",                 "HorseheadSkyViewFixed.wl"},
  {"horsehead_alternative.py",                   "HorseheadAlternative.wl"},
  {"horsehead_complete_demo.py",                 "HorseheadCompleteDemo.wl"},
  {"astronomy_viewer_gui.py",                    "AstronomyViewerGUI.wl"},
  {"interactive_sky_map.py",                     "InteractiveSkyMap.wl"}},
  Frame -> All,
  Background -> {None, {{LightGray}, White}},
  Alignment -> {Left, Center},
  Spacings -> {1.2, 0.6}];

cells = {

 Cell["From Python to Wolfram: Exploring NASA SkyView", "Title"],

 Cell["A walkthrough of the AstronomyImages toolkit port",
      "Subtitle"],

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
   "The open-source project ",
   StyleBox["AstronomyImages", "Input"],
   " is a small Python toolkit for pulling telescope imagery from the \
NASA SkyView Virtual Observatory, analysing the returned FITS files, \
and exposing the whole workflow through a couple of Tkinter GUIs. \
This notebook is a companion to the complete Wolfram-Language port \
shipped in the ",
   StyleBox["Wolfram/", "Input"],
   " folder of the same repository. It shows how the Python data \
pipeline maps to Wolfram idioms (",
   StyleBox["URLRead", "Input"], ", ",
   StyleBox["Import[_, \"FITS\"]", "Input"], ", ",
   StyleBox["ArrayPlot", "Input"], ", ",
   StyleBox["Manipulate", "Input"],
   "), downloads a real 0.5\[Degree] cutout of the Horsehead Nebula \
from NASA SkyView live, and lets you explore three surveys \
\[LongDash] DSS2 Red, 2MASS-J, and a physics-based simulation \
\[LongDash] in an interactive panel. Everything is public-domain \
data; the code is MIT-licensed."}],
      "Text"],

 imgCell["horsehead_DSS2_Red.png",
  "Figure 1. 600\[Times]600 DSS2 Red cutout of the Horsehead Nebula \
(Barnard 33) centred on RA 85.25\[Degree], Dec \[Minus]2.46\[Degree], \
0.5\[Degree] on a side. The image was fetched live from NASA \
SkyView and rendered by ArrayPlot. The bright star at the upper \
left is \[Sigma] Orionis; the dark silhouette is the nebula itself, \
a foreground dust cloud absorbing the emission of IC 434 behind it."],

 (* ---- Section 1 ---- *)

 Cell["1. What the original toolkit does", "Section"],

 Cell[TextData[{"The Python project in the parent directory contains \
five kinds of scripts:\n",
   "\[Bullet] downloading FITS cutouts from the ",
   ButtonBox["NASA SkyView Virtual Observatory",
     BaseStyle -> "Hyperlink",
     ButtonData -> {URL["https://skyview.gsfc.nasa.gov/"], None}],
   " (public domain);\n",
   "\[Bullet] fallback JPEG queries against SDSS SkyServer and STScI DSS2;\n",
   "\[Bullet] a physics-based simulation of the Horsehead dark nebula;\n",
   "\[Bullet] two Tkinter GUIs \[LongDash] a simple viewer and a \
clickable all-sky map;\n",
   "\[Bullet] analysis scaffolding (histograms, percentile clipping, \
optical-depth estimates).\n\n",
   "The Wolfram port lives in ",
   StyleBox["Wolfram/", "Input"],
   " and mirrors every script as a self-contained ",
   StyleBox[".wl", "Input"], " package."}],
      "Text"],

 Cell[BoxData[ToBoxes[cheatGrid]],
  "Output", ShowCellLabel -> False,
  CellMargins -> {{0, 0}, {8, 14}}],

 Cell["1.1 Python\[LongDash]to\[LongDash]Wolfram cheat-sheet",
      "Subsection"],

 Cell[TextData[{"A few non-obvious translations:\n",
   "\[Bullet] ",
   StyleBox["requests.get(url, params=..)", "Input"],
   " \[Rule] ",
   StyleBox["URLRead[HTTPRequest[url<>\"?\"<>URLQueryEncode[params]]][\"BodyByteArray\"]",
     "Input"], "\n",
   "\[Bullet] ",
   StyleBox["astropy.io.fits.open(..)", "Input"],
   " \[Rule] ",
   StyleBox["Import[file, {\"FITS\", \"Data\"}]", "Input"],
   " (+ ",
   StyleBox["\"Metadata\"", "Input"], ")\n",
   "\[Bullet] ",
   StyleBox["numpy.linspace(a,b,n)", "Input"],
   " \[Rule] ",
   StyleBox["Subdivide[a, b, n-1]", "Input"], "\n",
   "\[Bullet] ",
   StyleBox["numpy.percentile(a, q)", "Input"],
   " \[Rule] ",
   StyleBox["Quantile[a, q/100]", "Input"], "\n",
   "\[Bullet] ",
   StyleBox["numpy.random.normal / poisson", "Input"],
   " \[Rule] ",
   StyleBox["RandomVariate[NormalDistribution[..]]", "Input"],
   " / ",
   StyleBox["RandomVariate[PoissonDistribution[..]]", "Input"], "\n",
   "\[Bullet] ",
   StyleBox["matplotlib.imshow + colorbar", "Input"],
   " \[Rule] ",
   StyleBox["ArrayPlot[Reverse[data], PlotLegends -> Automatic]",
     "Input"], "\n",
   "\[Bullet] ",
   StyleBox["Tkinter widget tree", "Input"],
   " \[Rule] ",
   StyleBox["DynamicModule[..] + Manipulate[..]", "Input"], "\n",
   "\[Bullet] ",
   StyleBox["json.dump", "Input"],
   " \[Rule] ",
   StyleBox["Export[file, assoc, \"JSON\"]", "Input"]}],
      "Text"],

 Cell[TextData[{"A subtlety worth flagging: ",
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
working idiom is ",
   StyleBox["URLRead[HTTPRequest[url]][\"BodyByteArray\"]", "Input"],
   ", which returns the raw bytes."}],
      "Text"],

 (* ---- Section 2 ---- *)

 Cell["2. Downloading a real image from NASA SkyView", "Section"],

 Cell[TextData[{"Everything SkyView needs is packed into the query \
string. The cell below fetches a 0.5\[Degree] cutout around the \
Horsehead Nebula (Barnard 33) in the DSS2 Red band and writes the \
raw FITS file next to the notebook."}],
      "Text"],

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
skyViewBytes[ra_, dec_, survey_, fov_:0.5, pixels_:600] :=\n\
  Module[{resp = URLRead[HTTPRequest[skyViewURL[ra, dec, survey, fov, pixels]]]},\n\
    If[resp === $Failed || resp[\"StatusCode\"] =!= 200,\n\
      $Failed, resp[\"BodyByteArray\"]]];"],
      "Input"],

 Cell[TextData[{"The notebook ships with three pre-downloaded cutouts \
in ",
   StyleBox["images/", "Input"],
   ", so every cell below can be re-run offline. The cell below \
loads the DSS2 Red cutout and turns it into a numeric matrix. The \
FITS importer returns a one-entry Association \
(",
   StyleBox["<|1 -> array|>", "Input"],
   ") for single-HDU files and a longer Association for multi-HDU \
files; ",
   StyleBox["First[Values[..]]", "Input"],
   " picks the primary HDU either way."}], "Text"],

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
   " gives a publication-quality image in a single line."}], "Text"],

 Cell[BoxData["\
quickView[d_, cf_:\"GrayTones\"] :=\n\
  Module[{flat = Flatten[d], lo, hi},\n\
    {lo, hi} = Quantile[flat, {0.01, 0.995}];\n\
    ArrayPlot[Reverse[d],\n\
      ColorFunction -> ColorData[cf],\n\
      PlotRange -> {lo, hi},\n\
      ImageSize -> 480, Frame -> False]];\n\
\n\
quickView[data]"], "Input"],

 imgCell["horsehead_2MASS_J.png",
  "Figure 2. The same region in 2MASS-J (1.25 \[Mu]m, near-infrared). \
Infrared light penetrates dust more readily than visible light, so the \
Horsehead silhouette softens and background stars that were hidden \
in DSS2 Red punch through."],

 (* ---- Section 3 ---- *)

 Cell["3. Rudimentary image analysis", "Section"],

 Cell[TextData[{"This section shows the kind of lightweight analysis \
you can do on an imported FITS image with built-in Wolfram functions \
\[LongDash] no specialised toolbox required."}], "Text"],

 Cell["3.1 Summary statistics", "Subsection"],

 Cell[BoxData["flat = Flatten[data];"], "Input"],

 Cell[BoxData["\
<|\"Min\"    -> Min[flat],\n\
  \"Max\"    -> Max[flat],\n\
  \"Mean\"   -> N @ Mean[flat],\n\
  \"Median\" -> N @ Median[flat],\n\
  \"StdDev\" -> N @ StandardDeviation[flat]|>"], "Input"],

 Cell["3.2 Intensity histogram", "Subsection"],

 Cell[BoxData["\
Histogram[flat, 80,\n\
  ChartStyle -> RGBColor[0.27, 0.51, 0.71],\n\
  Frame -> True,\n\
  FrameLabel -> {\"Pixel intensity (ADU)\", \"Count\"},\n\
  PlotLabel -> \"Intensity distribution\",\n\
  ImageSize -> 520]"], "Input"],

 Cell["3.3 Dark-nebula segmentation", "Subsection"],

 Cell[TextData[{"The Horsehead is a ",
   StyleBox["dark", FontSlant -> Italic],
   " nebula \[LongDash] a foreground dust cloud silhouetted against \
the brighter H\[Alpha] emission region IC 434. A crude but effective \
segmentation is to flag every pixel below the 30th percentile."}],
      "Text"],

 Cell[BoxData["threshold = Quantile[flat, 0.30]"], "Input"],

 Cell[BoxData["darkMask = UnitStep[threshold - data];"], "Input"],

 Cell[BoxData["\
{quickView[data],\n\
 ArrayPlot[Reverse[darkMask],\n\
   ColorRules -> {0 -> White, 1 -> Red},\n\
   ImageSize -> 480, Frame -> False]}"], "Input"],

 Cell[TextData[{"Estimating the mean optical depth across the dark \
region gives a physically meaningful number. The Horsehead is ~0.5 \
mag optically thick in DSS2 Red, so we expect ",
   Cell[BoxData[FormBox[RowBox[{"\[Tau]", "\[TildeTilde]",
     "0.5"}], TraditionalForm]]],
   "."}], "Text"],

 Cell[BoxData["\
bgMean  = Mean[Pick[flat, Flatten[darkMask], 0]];\n\
nebMean = Mean[Pick[flat, Flatten[darkMask], 1]];\n\
tauEst  = -Log[nebMean / bgMean];\n\
<|\"bg (ADU)\"       -> N @ bgMean,\n\
  \"dark (ADU)\"     -> N @ nebMean,\n\
  \"optical depth\"  -> N @ tauEst,\n\
  \"extinction (mag)\" -> N[1.086 tauEst]|>"], "Input"],

 imgCell["horsehead_analysis_panel.png",
  "Figure 3. Four-panel summary of the DSS2 Red cutout: raw \
percentile-clipped image (upper left), log10-scaled version \
highlighting faint IC 434 emission (upper right), dark-nebula \
segmentation at the 30th percentile (lower left), and intensity \
histogram with the threshold marked (lower right)."],

 Cell["3.4 Edge detection, smoothing, log scaling", "Subsection"],

 Cell[TextData[{"Wolfram's image-processing functions operate \
directly on the 2-D array \[LongDash] no conversion needed."}], "Text"],

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

 (* ---- Section 4 ---- *)

 Cell["4. An interactive survey explorer", "Section"],

 Cell[TextData[{"The ",
   StyleBox["Manipulate", "Input"],
   " below turns the analysis into a live tool: scrub through three \
bundled surveys, change the colour map, slide the dark-threshold \
percentile up or down and watch the segmentation update in real \
time. No network round-trip is required \[LongDash] everything reads \
from the pre-downloaded files in ",
   StyleBox["images/", "Input"], "."}],
      "Text"],

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
  If[StringEndsQ[f, \".mx\"], Import[f],\n\
     loadFITSArray[f]];"], "Input"],

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

 imgCell["horsehead_simulated.png",
  "Figure 4. The simulated Horsehead included in the bundle: a \
diffuse background modulated by a low-frequency Gaussian blob, \
overlaid with two elongated Gaussians playing the role of the dust \
column. The pixel values are Poisson-sampled to mimic photon \
statistics. Having an offline, network-free image in the same \
framework makes the notebook run in any environment."],

 Cell[TextData[{"For the full all-sky clickable viewer, open ",
   StyleBox["Wolfram/InteractiveSkyMap.wl", "Input"],
   " in a notebook and evaluate ",
   StyleBox["InteractiveSkyMap`launchInteractiveSkyMap[]", "Input"],
   "."}], "Text"],

 (* ---- Section 5 ---- *)

 Cell["5. Takeaways", "Section"],

 Cell[TextData[{
   "\[Bullet] ",
   StyleBox["URLRead + HTTPRequest + URLQueryEncode", "Input"],
   " replace the Python ",
   StyleBox["requests", "Input"],
   " / ",
   StyleBox["urllib", "Input"],
   " dance in two lines.\n",
   "\[Bullet] ",
   StyleBox["Import[_, {\"FITS\", \"Data\"}]", "Input"],
   " reads SkyView's raw bytes with no ",
   StyleBox["BytesIO", "Input"],
   " shim \[LongDash] the only thing to watch out for is that it \
returns a single-entry Association for a single-HDU file.\n",
   "\[Bullet] ",
   StyleBox["ArrayPlot", "Input"], ", ",
   StyleBox["Histogram", "Input"], ", ",
   StyleBox["GradientFilter", "Input"],
   " and ",
   StyleBox["GaussianFilter", "Input"],
   " cover most of matplotlib + scipy.ndimage for this use case.\n",
   "\[Bullet] ",
   StyleBox["DynamicModule", "Input"],
   " + ",
   StyleBox["EventHandler", "Input"],
   " + ",
   StyleBox["Manipulate", "Input"],
   " replace Tkinter with a notebook-native UI you can evaluate on \
any platform.\n",
   "\[Bullet] Images shipped with the notebook come from NASA \
SkyView (public domain) and a built-in physics-based simulation."}],
      "Text"],

 Cell["6. References & licensing", "Section"],

 Cell[TextData[{
   "\[Bullet] ",
   ButtonBox["NASA SkyView Virtual Observatory",
     BaseStyle -> "Hyperlink",
     ButtonData -> {URL["https://skyview.gsfc.nasa.gov/"], None}],
   " \[LongDash] public domain, unrestricted use.\n",
   "\[Bullet] ",
   ButtonBox["Digitized Sky Survey (STScI)",
     BaseStyle -> "Hyperlink",
     ButtonData -> {URL["https://archive.stsci.edu/dss/acknowledging.html"],
                    None}],
   " \[LongDash] free for research and education with attribution.\n",
   "\[Bullet] ",
   ButtonBox["2MASS (IPAC/Caltech)",
     BaseStyle -> "Hyperlink",
     ButtonData -> {URL["https://www.ipac.caltech.edu/2mass/overview/about2mass.html"],
                    None}],
   " \[LongDash] public domain.\n",
   "\[Bullet] All code in this post is MIT-licensed.\n",
   "\[Bullet] Repository: ",
   ButtonBox["github.com/mthiel74/AstronomyImages",
     BaseStyle -> "Hyperlink",
     ButtonData -> {URL["https://github.com/mthiel74/AstronomyImages"],
                    None}], "."}],
      "Text"]
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
