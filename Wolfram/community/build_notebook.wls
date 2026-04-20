#!/usr/bin/env wolframscript
(* Build Wolfram/community/AstronomyImagesPost.nb \[LongDash] a Wolfram
   Community post on pulling real telescope imagery from NASA SkyView
   with the Wolfram Language. Ships a dozen cutouts of famous targets,
   a 2MASS J/H/K false-colour composite, a multi-wavelength Horsehead
   panel, and a live downloader you can aim at any RA/Dec. *)

hereDir = DirectoryName[If[$ScriptCommandLine =!= {}, $InputFileName,
                           NotebookFileName[]]];
imgDir  = FileNameJoin[{hereDir, "images"}];

imgWidth = 520;

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
                   imgWidth - 8, Alignment -> Center]},
             Alignment -> Center, Spacings -> 0.6]]]],
  "Output", ShowCellLabel -> False,
  TextAlignment -> Center,
  CellMargins -> {{0, 0}, {10, 14}}];

cells = {

 Cell["Exploring NASA SkyView with the Wolfram Language", "Title"],

 Cell["Live telescope imagery, dark-nebula analysis, and a \
catalog-driven sky-survey downloader", "Subtitle"],

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
The friction is entirely in the tooling."}], "Text"],

 Cell[TextData[{"This notebook shows how little Wolfram Language \
code that friction actually needs. One ",
   StyleBox["URLRead", "Input"], " fetches a FITS cutout; one ",
   StyleBox["Import", "Input"], " decodes it; one ",
   StyleBox["ArrayPlot", "Input"], " renders it; one ",
   StyleBox["DynamicModule", "Input"],
   " wraps the lot into a live multi-survey explorer. Along the way \
we pull a 0.5\[Degree] DSS2 Red cutout of the Horsehead Nebula, build \
a 2MASS J/H/K false-colour composite with an asinh stretch, tour six \
Messier targets (M1, M16, M31, M42, M45, M51) across multiple colour \
maps, and expose a catalog-driven downloader that lets you aim the \
pipeline at any RA/Dec on the sky."}], "Text"],

 Cell[TextData[{"The notebook is self-contained: a single \
initialisation cell at the top downloads the FITS data it needs \
directly from NASA SkyView on first run, caches it next to the \
notebook, and every subsequent open is offline. No other files to \
attach."}], "Text"],

 Cell[TextData[{"The code also ships as ",
   StyleBox[".wl", "Input"],
   " packages in the ",
   StyleBox["Wolfram/", "Input"],
   " folder of the ",
   ButtonBox["AstronomyImages",
     BaseStyle -> "Hyperlink",
     ButtonData -> {URL["https://github.com/mthiel74/AstronomyImages"],
                    None}],
   " repository (MIT-licensed); a parallel Python implementation \
lives at the repository root for readers who prefer that stack."}],
      "Text"],

 imgCell["horsehead_DSS2_Red.png",
  "Figure 1. A live 600\[Times]600 DSS2 Red cutout of the Horsehead \
Nebula (Barnard 33), centred at RA 85.25\[Degree], \
Dec \[Minus]2.46\[Degree], 0.5\[Degree] on a side. The bright star at \
the upper left is \[Sigma] Orionis; the horse-head silhouette is a \
foreground molecular cloud eating into the emission nebula IC 434 \
behind it. Rendered by Wolfram in a single ArrayPlot call. \
DSS2 / STScI."],

 (* ---- Initialisation (self-contained bootstrapping) ---- *)

 Cell["Initialisation", "Section"],

 Cell[TextData[{"The cell below is the only thing you have to \
evaluate before the rest of the post runs. It defines the shared \
helpers (",
   StyleBox["skyViewURL", "Input"], ", ",
   StyleBox["skyViewBytes", "Input"], ", ",
   StyleBox["loadFITSArray", "Input"],
   ") and, the first time you evaluate it, downloads a bundle of \
FITS cutouts (~18 MB total) from NASA SkyView into an ",
   StyleBox["images/", "Input"],
   " folder next to this notebook. Subsequent opens reuse the \
cached files, so the notebook then runs fully offline. The cell \
is flagged as an ",
   StyleBox["InitializationCell", "Input"],
   ", so Wolfram will offer to evaluate it automatically when the \
notebook is first opened."}], "Text"],

 Cell[BoxData["\
imagesDir = FileNameJoin[{NotebookDirectory[], \"images\"}];\n\
If[!DirectoryQ[imagesDir], CreateDirectory[imagesDir]];\n\
\n\
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
      $Failed, resp[\"BodyByteArray\"]]];\n\
\n\
loadFITSArray[file_] :=\n\
  With[{raw = Import[file, {\"FITS\", \"Data\"}]},\n\
    Which[\n\
      AssociationQ[raw],                First[Values[raw]],\n\
      ListQ[raw] && Depth[raw] >= 4,    First[raw],\n\
      True,                             raw]];\n\
\n\
cachedFITS[fileName_, ra_, dec_, survey_, fov_: 0.5, pixels_: 600] :=\n\
  Module[{path = FileNameJoin[{imagesDir, fileName}], bytes},\n\
    If[FileExistsQ[path] && FileByteCount[path] > 8000, Return[path]];\n\
    Print[\"downloading \", fileName, \" (\", survey, \", \", fov, \"\[Degree])\[Ellipsis]\"];\n\
    bytes = skyViewBytes[ra, dec, survey, fov, pixels];\n\
    If[!ByteArrayQ[bytes] || Length[bytes] < 8000,\n\
      Print[\"  ! failed\"]; Return[$Failed]];\n\
    BinaryWrite[path, bytes]; Close[path];\n\
    path];\n\
\n\
(* Horsehead: five bands, one position *)\n\
cachedFITS[\"horsehead_DSS2_Red.fits\", 85.2479, -2.4583, \"DSS2 Red\", 0.5, 600];\n\
cachedFITS[\"horsehead_2MASS_J.fits\",  85.2479, -2.4583, \"2MASS-J\",  0.5, 600];\n\
cachedFITS[\"horsehead_2MASS_H.fits\",  85.2479, -2.4583, \"2MASS-H\",  0.5, 600];\n\
cachedFITS[\"horsehead_2MASS_K.fits\",  85.2479, -2.4583, \"2MASS-K\",  0.5, 600];\n\
cachedFITS[\"horsehead_WISE_3p4.fits\", 85.2479, -2.4583, \"WISE 3.4\", 0.5, 600];\n\
\n\
(* Six Messier targets in DSS2 Red *)\n\
cachedFITS[\"M42_Orion_Nebula_DSS2_Red.fits\",  83.8213, -5.3903, \"DSS2 Red\", 1.0,  600];\n\
cachedFITS[\"M45_Pleiades_DSS2_Red.fits\",     56.870,   24.105, \"DSS2 Red\", 2.5,  600];\n\
cachedFITS[\"M31_Andromeda_DSS2_Red.fits\",    10.685,   41.269, \"DSS2 Red\", 3.0,  600];\n\
cachedFITS[\"M16_Eagle_Nebula_DSS2_Red.fits\", 274.700, -13.807, \"DSS2 Red\", 0.5,  600];\n\
cachedFITS[\"M51_Whirlpool_DSS2_Red.fits\",    202.469,  47.196, \"DSS2 Red\", 0.2,  600];\n\
cachedFITS[\"M1_Crab_Nebula_DSS2_Red.fits\",   83.633,   22.014, \"DSS2 Red\", 0.15, 600];\n\
\n\
(* Simulated Horsehead (pure Wolfram, no network) *)\n\
With[{simPath = FileNameJoin[{imagesDir, \"horsehead_simulated.mx\"}]},\n\
  If[!FileExistsQ[simPath],\n\
    Print[\"generating simulated Horsehead\[Ellipsis]\"];\n\
    SeedRandom[42];\n\
    Module[{size = 600, xs, X, Y, bg, mask, tau, img, sim},\n\
      xs = Subdivide[-1., 1., size - 1];\n\
      X = Table[x, {y, xs}, {x, xs}];\n\
      Y = Table[y, {y, xs}, {x, xs}];\n\
      bg = 100. + 150. Exp[-(X^2 + Y^2)/0.5];\n\
      Do[bg = bg + RandomReal[{10, 30}]\n\
                  Sin[RandomReal[{5, 15}] (X Cos[RandomReal[{0, Pi}]]\n\
                                        + Y Sin[RandomReal[{0, Pi}]])], {5}];\n\
      mask = Exp[-(((X - 0.10)/0.30)^2 + ((Y + 0.20)/0.40)^2) 8.]\n\
           + 0.7 Exp[-(((X - 0.15)/0.15)^2 + ((Y - 0.15)/0.35)^2) 6.];\n\
      tau = 3.0 mask;\n\
      img = bg Exp[-tau];\n\
      img = Map[RandomVariate[PoissonDistribution[Max[#, 0.]]] &, img, {2}];\n\
      sim = Clip[img + RandomVariate[NormalDistribution[0.,\n\
              0.02 Mean[Flatten[img]]], {size, size}], {0., Infinity}];\n\
      Export[simPath, sim]]]];\n\
\n\
Print[\"Initialisation complete. Bundled files:\"];\n\
FileNames[\"*.fits\" | \"*.mx\", imagesDir]"], "Input",
  InitializationCell -> True],

 (* ---- Section 1: Downloading ---- *)

 Cell["1. Downloading a real image from NASA SkyView", "Section"],

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

 Cell[TextData[{StyleBox["Heads-up on URLRead vs URLExecute.",
     FontWeight -> Bold],
   " ",
   StyleBox["URLExecute[url, \"ByteArray\"]", "Input"],
   " looks like it should give you the raw response bytes \[LongDash] \
it doesn't. The second argument of ",
   StyleBox["URLExecute", "Input"],
   " is the ",
   StyleBox["Import", "Input"],
   " element name, not the response encoding. When SkyView returns ",
   StyleBox["Content-Type: image/fits", "Input"],
   " the call therefore tries to parse the body as a FITS file and \
extract a field called \"ByteArray\", which does not exist. The \
working idiom is ",
   StyleBox["URLRead[HTTPRequest[url]][\"BodyByteArray\"]", "Input"],
   " as above."}], "Text"],

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

 Cell["2. One object, four wavelengths", "Section"],

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

 Cell[TextData[{"The cell below loads four bands from the bundled \
FITS files and builds the 2\[Times]2 panel in one Grid. Every panel \
is a single ",
   StyleBox["ArrayPlot", "Input"],
   " with a per-band colour map and percentile-clipped contrast."}],
      "Text"],

 Cell[BoxData["\
waveBands = {\n\
  {\"DSS2 Red (656 nm)\",    \"horsehead_DSS2_Red.fits\",  \"GrayTones\"},\n\
  {\"2MASS J (1.25 \[Mu]m)\",    \"horsehead_2MASS_J.fits\",   \"GrayTones\"},\n\
  {\"2MASS K (2.17 \[Mu]m)\",    \"horsehead_2MASS_K.fits\",   \"SunsetColors\"},\n\
  {\"WISE W1 (3.4 \[Mu]m)\",     \"horsehead_WISE_3p4.fits\",  \"SunsetColors\"}};\n\
\n\
wavePanels = Table[\n\
  With[{d = loadFITSArray[FileNameJoin[\n\
              {NotebookDirectory[], \"images\", b[[2]]}]]},\n\
    Labeled[\n\
      ArrayPlot[Reverse[d],\n\
        ColorFunction -> ColorData[b[[3]]],\n\
        PlotRange -> Quantile[Flatten[d], {0.01, 0.995}],\n\
        Frame -> False, ImageSize -> 250],\n\
      Style[b[[1]], Bold, 12]]],\n\
  {b, waveBands}];\n\
\n\
Labeled[Grid[Partition[wavePanels, 2], Spacings -> 1],\n\
  Style[\"Horsehead Nebula: one object, four wavelengths\",\n\
        Bold, 14], Top]"], "Input"],

 Cell[TextData[{"The pre-rendered PNG below is the exact output of \
that cell, embedded so readers on GitHub or the Community post see \
the figure without evaluating anything."}], "Text"],

 imgCell["horsehead_multiwavelength.png",
  "Figure 2. The Horsehead Nebula in four bands spanning a factor of \
~5000 in wavelength. DSS2 Red (656 nm) shows the classic silhouette \
absorbing H\[Alpha]. 2MASS J (1.25 \[Mu]m) is only 2\[Times] redward \
but the extinction is already half its optical value \[LongDash] \
field stars punch through. At 2MASS K (2.17 \[Mu]m) the nebula is \
almost transparent. WISE W1 (3.4 \[Mu]m, shown in SunsetColors) \
reveals warm dust emission from the nebula itself, plus the \
broadband glow of PAH molecules along the IC 434 ridge. \
DSS2 / STScI; 2MASS / IPAC; WISE / NASA."],

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

 Cell["3. A 2MASS J/H/K false-colour composite", "Section"],

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
while preserving low-surface-brightness structure. 2MASS / IPAC."],

 (* ---- Section 6: rudimentary analysis ---- *)

 Cell["4. Rudimentary image analysis", "Section"],

 Cell[TextData[{"Every tool you need for a first-pass numerical look \
at a FITS image is in the top-level ",
   StyleBox["System`", "Input"],
   " context. No toolbox imports."}], "Text"],

 Cell["4.1 Summary statistics", "Subsection"],

 Cell[BoxData["flat = Flatten[data];"], "Input"],

 Cell[BoxData["\
<|\"Min\"    -> Min[flat],\n\
  \"Max\"    -> Max[flat],\n\
  \"Mean\"   -> N @ Mean[flat],\n\
  \"Median\" -> N @ Median[flat],\n\
  \"StdDev\" -> N @ StandardDeviation[flat]|>"], "Input"],

 Cell["4.2 Intensity histogram", "Subsection"],

 Cell[BoxData["\
Histogram[flat, 80,\n\
  ChartStyle -> RGBColor[0.27, 0.51, 0.71],\n\
  Frame -> True,\n\
  FrameLabel -> {\"Pixel intensity (ADU)\", \"Count\"},\n\
  PlotLabel -> \"Intensity distribution\",\n\
  ImageSize -> 520]"], "Input"],

 Cell["4.3 Dark-nebula segmentation and optical depth", "Subsection"],

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
threshold marked (lower right). DSS2 / STScI."],

 Cell["4.4 Smoothing, log scaling, edge detection", "Subsection"],

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

 Cell["5. A Messier gallery", "Section"],

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

 Cell[TextData[{"The cell below holds the gallery as a catalog, \
loads each bundled FITS, and arranges the six panels in a 3\[Times]2 \
Grid. Each target gets its own ColorFunction \[LongDash] \
",
   StyleBox["Rainbow", "Input"], " for M42, ",
   StyleBox["DeepSeaColors", "Input"], " for M45, and so on \
\[LongDash] so the gallery doubles as a palette sampler."}], "Text"],

 Cell[BoxData["\
messierGallery = <|\n\
  \"M42 (Orion Nebula)\" -> <|\"file\" -> \"M42_Orion_Nebula_DSS2_Red.fits\",\n\
                             \"cmap\" -> \"Rainbow\",\n\
                             \"caption\" -> \"Stellar nursery, 1344 ly away\"|>,\n\
  \"M45 (Pleiades)\"     -> <|\"file\" -> \"M45_Pleiades_DSS2_Red.fits\",\n\
                             \"cmap\" -> \"DeepSeaColors\",\n\
                             \"caption\" -> \"Open cluster, ~110 Myr old\"|>,\n\
  \"M31 (Andromeda)\"    -> <|\"file\" -> \"M31_Andromeda_DSS2_Red.fits\",\n\
                             \"cmap\" -> \"SunsetColors\",\n\
                             \"caption\" -> \"Nearest spiral galaxy, 2.5 Mly\"|>,\n\
  \"M16 (Eagle Nebula)\" -> <|\"file\" -> \"M16_Eagle_Nebula_DSS2_Red.fits\",\n\
                             \"cmap\" -> \"TemperatureMap\",\n\
                             \"caption\" -> \"Pillars of Creation, HST target\"|>,\n\
  \"M51 (Whirlpool)\"    -> <|\"file\" -> \"M51_Whirlpool_DSS2_Red.fits\",\n\
                             \"cmap\" -> \"AuroraColors\",\n\
                             \"caption\" -> \"Grand-design spiral, 23 Mly\"|>,\n\
  \"M1 (Crab Nebula)\"   -> <|\"file\" -> \"M1_Crab_Nebula_DSS2_Red.fits\",\n\
                             \"cmap\" -> \"SolarColors\",\n\
                             \"caption\" -> \"Supernova remnant of 1054 AD\"|>|>;\n\
\n\
galleryPanels = Table[\n\
  With[{meta = messierGallery[k],\n\
        d = loadFITSArray[FileNameJoin[\n\
              {NotebookDirectory[], \"images\", messierGallery[k][\"file\"]}]]},\n\
    Labeled[\n\
      ArrayPlot[Reverse[d],\n\
        ColorFunction -> ColorData[meta[\"cmap\"]],\n\
        PlotRange -> Quantile[Flatten[d], {0.01, 0.995}],\n\
        Frame -> False, ImageSize -> 220],\n\
      Column[{Style[k, Bold, 11],\n\
              Style[meta[\"caption\"], Italic, FontSize -> 9, Gray]},\n\
             Alignment -> Center]]],\n\
  {k, Keys[messierGallery]}];\n\
\n\
Grid[Partition[galleryPanels, 3], Spacings -> {1, 1.5}]"], "Input"],

 Cell[TextData[{"The PNG below is the same computation baked to \
disk by ",
   StyleBox["make_assets.wls", "Input"],
   " so readers of the static PDF / GitHub render see the figure \
without having to evaluate anything."}], "Text"],

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
DSS2 / STScI."],

 Cell[TextData[{"For individual high-resolution plates of each \
target, see ",
   StyleBox["images/gallery_*.png", "Input"],
   " in the notebook directory, or re-run ",
   StyleBox["make_assets.wls", "Input"],
   " to regenerate them from scratch."}], "Text"],

 (* ---- Section 8: interactive explorer ---- *)

 Cell["6. Interactive explorers", "Section"],

 Cell["6.1 Offline Manipulate over the bundled surveys", "Subsection"],

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

 Cell["6.2 Live downloader \[LongDash] any object, any survey",
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
reports success, image size, and bytes received. DSS2 / STScI."],

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

 Cell["6.3 Click-anywhere sky map", "Subsection"],

 Cell[TextData[{"Clicking a Wolfram ",
   StyleBox["Graphics", "Input"],
   " object and getting the coordinates out is what ",
   StyleBox["EventHandler", "Input"],
   " + ",
   StyleBox["MousePosition[\"Graphics\"]", "Input"],
   " are for. Wrap an all-sky Mercator in that handler, stash the \
click in ",
   StyleBox["{selRA, selDec}", "Input"],
   ", and pipe those straight into ",
   StyleBox["skyViewBytes", "Input"],
   " when the Download button is pressed. The cell below is about \
forty lines end-to-end."}], "Text"],

 Cell[BoxData["\
skyLandmarks = <|\n\
  \"Horsehead\"     -> {85.25,  -2.46},\n\
  \"M42 Orion\"     -> {83.82,  -5.39},\n\
  \"M45 Pleiades\"  -> {56.87,  24.11},\n\
  \"M31 Andromeda\" -> {10.68,  41.27},\n\
  \"M16 Eagle\"     -> {274.70, -13.81},\n\
  \"M51 Whirlpool\" -> {202.47, 47.20},\n\
  \"M1 Crab\"       -> { 83.63, 22.01}|>;\n\
\n\
skyMap[selRA_, selDec_, fov_] := Graphics[{\n\
    (* Milky Way approximation *)\n\
    {Purple, Opacity[0.15], Thickness[0.03],\n\
      Line[{#, -29 + 60 Sin[(# - 266) Degree]} & /@ Subdivide[0., 360., 99]]},\n\
    (* Celestial equator + ecliptic *)\n\
    {Blue, Opacity[0.4], Line[{{0, 0}, {360, 0}}]},\n\
    {Orange, Opacity[0.4], Dashed,\n\
      Line[{#, 23.4 Sin[# Degree]} & /@ Subdivide[0., 360., 99]]},\n\
    (* Landmarks *)\n\
    KeyValueMap[\n\
      Function[{n, p},\n\
        {Red, PointSize[0.012], Point[p],\n\
         Text[Style[n, 9, Darker[Red]], p, {0, -1.5}]}],\n\
      skyLandmarks],\n\
    (* Selection \[LongDash] square footprint of side fov\[Degree] *)\n\
    If[NumericQ[selRA],\n\
      {Green, PointSize[0.02], Point[{selRA, selDec}],\n\
       (* True-size FOV square (vanishingly small for FOV < ~1\[Degree]) *)\n\
       Thickness[0.003],\n\
       Line[{{selRA - fov/2, selDec - fov/2},\n\
             {selRA + fov/2, selDec - fov/2},\n\
             {selRA + fov/2, selDec + fov/2},\n\
             {selRA - fov/2, selDec + fov/2},\n\
             {selRA - fov/2, selDec - fov/2}}],\n\
       (* Dashed \"minimum visible\" outline so slider feedback is \n\
          always obvious, even for sub-degree FOV *)\n\
       Dashed, Opacity[0.6],\n\
       With[{r = Max[fov/2, 4]},\n\
        Line[{{selRA - r, selDec - r}, {selRA + r, selDec - r},\n\
              {selRA + r, selDec + r}, {selRA - r, selDec + r},\n\
              {selRA - r, selDec - r}}]],\n\
       (* Live FOV label next to the cursor *)\n\
       Opacity[1], Black,\n\
       Text[Style[\"FOV=\" <> ToString[NumberForm[fov, {3, 2}]] <> \"\[Degree]\",\n\
                  Bold, 10, Background -> RGBColor[1, 1, 1, 0.7]],\n\
            {selRA, selDec + Max[fov/2, 4] + 4}]}, {}]\n\
  }, PlotRange -> {{360, 0}, {-90, 90}},\n\
     Frame -> True,\n\
     FrameLabel -> {\"Right Ascension (\[Degree])\", \"Declination (\[Degree])\"},\n\
     AspectRatio -> 1/2, GridLines -> Automatic,\n\
     GridLinesStyle -> Directive[Gray, Opacity[0.3], Dashed],\n\
     ImageSize -> 680];\n\
\n\
(* Helper: fetch + decode one band, return {bytes, array} or $Failed *)\n\
fetchBand[ra_, dec_, surv_, f_, px_] := Module[{b, tmp, a},\n\
  b = skyViewBytes[ra, dec, surv, f, px];\n\
  If[!ByteArrayQ[b] || Length[b] < 8000, Return[$Failed]];\n\
  tmp = CreateFile[]; BinaryWrite[tmp, b]; Close[tmp];\n\
  a = loadFITSArray[tmp]; DeleteFile[tmp];\n\
  If[MatrixQ[a, NumericQ], {b, a}, $Failed]];\n\
\n\
bandCmap = <|\"DSS2 Red\" -> \"GrayTones\", \"2MASS-J\" -> \"GrayTones\",\n\
             \"2MASS-K\" -> \"SunsetColors\", \"WISE 3.4\" -> \"SunsetColors\"|>;\n\
\n\
(* Build the display panel \[LongDash] shared by Dynamic render + Save PNG *)\n\
buildResultPanel[bd_Association, ra_, dec_, f_] :=\n\
  If[Length[bd] == 0,\n\
    Framed[Style[\"(no image downloaded yet)\", 12, Gray],\n\
      ImageSize -> {540, 300}, Alignment -> Center],\n\
    Module[{panels},\n\
      panels = KeyValueMap[\n\
        Function[{b, d},\n\
          With[{lohi = Quantile[Flatten[d], {0.01, 0.995}]},\n\
            Labeled[\n\
              ArrayPlot[Reverse[d],\n\
                ColorFunction -> ColorData[bandCmap[b]],\n\
                PlotRange -> lohi, Frame -> False,\n\
                ImageSize -> If[Length[bd] == 1, 540, 260]],\n\
              Style[b, Bold, 11]]]], bd];\n\
      If[Length[bd] == 1, First[panels],\n\
        Labeled[Grid[Partition[panels, 2, 2, 1, \"\"], Spacings -> 1],\n\
          Style[\"RA=\" <> ToString[NumberForm[ra, {5, 2}]] <>\n\
                \"\[Degree], Dec=\" <> ToString[NumberForm[dec, {5, 2}]] <>\n\
                \"\[Degree], FOV=\" <> ToString[NumberForm[f, {3, 2}]] <>\n\
                \"\[Degree]\", Bold, 12], Top]]];\n\
\n\
DynamicModule[{selRA = 85.25, selDec = -2.46, fov = 0.5,\n\
    survey = \"DSS2 Red\", allBands = False,\n\
    bandData = <||>, bandBytes = <||>,\n\
    status = \"Click the map to select coordinates.\"},\n\
  Column[{\n\
    EventHandler[Dynamic[skyMap[selRA, selDec, fov]],\n\
      {\"MouseClicked\" :> With[{p = MousePosition[\"Graphics\"]},\n\
        If[Head[p] === List && Length[p] == 2, {selRA, selDec} = p]]}],\n\
    Row[{\"RA: \",  Dynamic[NumberForm[selRA, {6, 3}]], \"\[Degree]   \",\n\
         \"Dec: \", Dynamic[NumberForm[selDec, {6, 3}]], \"\[Degree]   \",\n\
         \"FOV: \", Slider[Dynamic[fov], {0.05, 3, 0.05}],\n\
         \" \", Dynamic[NumberForm[fov, {3, 2}]], \"\[Degree]\"}],\n\
    Row[{\"Survey: \", PopupMenu[Dynamic[survey],\n\
           {\"DSS2 Red\", \"2MASS-J\", \"2MASS-K\", \"WISE 3.4\"}], \"   \",\n\
         Checkbox[Dynamic[allBands]], \" All spectral bands\"}],\n\
    Row[{\n\
      Button[Style[\"Download\", Bold, White],\n\
        bandData = <||>; bandBytes = <||>;\n\
        If[allBands,\n\
          Module[{bands = {\"DSS2 Red\", \"2MASS-J\", \"2MASS-K\", \"WISE 3.4\"},\n\
                  done = 0, r},\n\
            Do[status = \"Downloading \" <> b <> \" (\" <>\n\
                 ToString[done + 1] <> \"/\" <> ToString[Length[bands]] <> \")\[Ellipsis]\";\n\
               r = fetchBand[selRA, selDec, b, fov, 600];\n\
               If[r =!= $Failed,\n\
                 bandBytes[b] = r[[1]]; bandData[b] = r[[2]]];\n\
               done++, {b, bands}];\n\
            status = \"\[Checkmark] \" <> ToString[Length[bandData]] <> \"/\" <>\n\
                     ToString[Length[bands]] <> \" bands downloaded\"],\n\
          status = \"Contacting SkyView\[Ellipsis]\";\n\
          Module[{r = fetchBand[selRA, selDec, survey, fov, 600]},\n\
            If[r === $Failed,\n\
              status = \"\[Times] SkyView returned no data\",\n\
              bandBytes[survey] = r[[1]]; bandData[survey] = r[[2]];\n\
              status = \"\[Checkmark] \" <> ToString[Dimensions[r[[2]]]] <>\n\
                       \" from \" <> survey]]],\n\
        Background -> RGBColor[0.12, 0.45, 0.80], ImageSize -> 110],\n\
      Spacer[8],\n\
      Button[\"Save FITS\",\n\
        If[Length[bandBytes] == 0,\n\
          status = \"\[Times] Nothing to save \[LongDash] press Download first\",\n\
          Module[{p = SystemDialogInput[\"FileSave\", \"cutout.fits\"]},\n\
            If[StringQ[p],\n\
              If[Length[bandBytes] == 1,\n\
                BinaryWrite[p, First[Values[bandBytes]]]; Close[p];\n\
                status = \"Saved \" <> p,\n\
                Module[{dir = DirectoryName[p], base = FileBaseName[p],\n\
                        written = 0},\n\
                  KeyValueMap[Function[{b, bytes},\n\
                    Module[{q = FileNameJoin[{dir, base <> \"_\" <>\n\
                            StringReplace[b, {\" \" -> \"_\", \"-\" -> \"_\"}] <>\n\
                            \".fits\"}]},\n\
                      BinaryWrite[q, bytes]; Close[q]; written++]],\n\
                    bandBytes];\n\
                  status = \"Saved \" <> ToString[written] <>\n\
                           \" FITS files to \" <> dir]]]]],\n\
        ImageSize -> 100],\n\
      Spacer[4],\n\
      Button[\"Save PNG\",\n\
        If[Length[bandData] == 0,\n\
          status = \"\[Times] Nothing to save \[LongDash] press Download first\",\n\
          Module[{p = SystemDialogInput[\"FileSave\", \"cutout.png\"]},\n\
            If[StringQ[p],\n\
              Export[p, Rasterize[\n\
                buildResultPanel[bandData, selRA, selDec, fov],\n\
                ImageResolution -> 150]];\n\
              status = \"Saved \" <> p]]],\n\
        ImageSize -> 100]\n\
    }],\n\
    Style[Dynamic[status], Italic, FontSize -> 11],\n\
    Dynamic[buildResultPanel[bandData, selRA, selDec, fov]]\n\
  }, Spacings -> 1]]"], "Input"],

 imgCell["click_map_preview.png",
  "Figure 8. Static preview of \[Section]6.3 after a click near \
the Horsehead and a press of Download. The upper panel is the \
clickable all-sky Mercator (celestial equator in blue, ecliptic \
dashed in orange, Milky Way in purple, seven landmarks in red; \
the green marker is the cursor's current selection). The lower \
panel is the resulting DSS2 Red cutout. DSS2 / STScI."],

 Cell[TextData[{"A more elaborate version with constellation \
patches and a save-FITS/PNG button lives in ",
   StyleBox["Wolfram/InteractiveSkyMap.wl", "Input"],
   " (",
   StyleBox["InteractiveSkyMap`launchInteractiveSkyMap[]", "Input"],
   "), but the cell above already covers the main pedagogical \
point: ",
   StyleBox["EventHandler", "Input"],
   " + ",
   StyleBox["MousePosition[\"Graphics\"]", "Input"],
   " turns any 2-D plot into a coordinate picker."}], "Text"],

 (* ---- Section 7: simulated ---- *)

 Cell["7. When the network is down: a physics-based simulation",
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

 (* ---- Section 10: acknowledgements ---- *)

 Cell["8. Acknowledgements", "Section"],

 Cell[TextData[{"All images in this post were fetched from the ",
   ButtonBox["NASA SkyView Virtual Observatory",
     BaseStyle -> "Hyperlink",
     ButtonData -> {URL["https://skyview.gsfc.nasa.gov/"], None}],
   ". The standard acknowledgement sentences for the two surveys \
used are reproduced below."}], "Text"],

 Cell[TextData[{
   "\[OpenCurlyDoubleQuote]The Digitized Sky Surveys were produced at \
the Space Telescope Science Institute under U.S. Government grant \
NAG W-2166. The images of these surveys are based on photographic \
data obtained using the Oschin Schmidt Telescope on Palomar Mountain \
and the UK Schmidt Telescope.\[CloseCurlyDoubleQuote]"}], "Text"],

 Cell[TextData[{
   "\[OpenCurlyDoubleQuote]This publication makes use of data \
products from the Two Micron All Sky Survey, which is a joint project \
of the University of Massachusetts and the Infrared Processing and \
Analysis Center / California Institute of Technology, funded by the \
National Aeronautics and Space Administration and the National \
Science Foundation.\[CloseCurlyDoubleQuote]"}], "Text"],

 Cell[TextData[{"WISE data courtesy of NASA / JPL-Caltech / UCLA. \
Code in this repository is MIT-licensed."}], "Text"],

 (* ---- Section 11: takeaways ---- *)

 Cell["9. Takeaways", "Section"],

 Cell[TextData[{
   "\[Bullet] ",
   StyleBox["URLRead", "Input"], " + ",
   StyleBox["HTTPRequest", "Input"], " + ",
   StyleBox["URLQueryEncode", "Input"],
   " pull a SkyView FITS cutout in two lines.\n",
   "\[Bullet] ",
   StyleBox["Import[_, {\"FITS\", \"Data\"}]", "Input"],
   " decodes the bytes directly; for a single-HDU file it returns a \
one-entry Association, so ",
   StyleBox["First[Values[..]]", "Input"],
   " is the robust unwrap.\n",
   "\[Bullet] ",
   StyleBox["ArrayPlot", "Input"], ", ",
   StyleBox["Histogram", "Input"], ", ",
   StyleBox["GradientFilter", "Input"],
   " and ",
   StyleBox["GaussianFilter", "Input"],
   " operate directly on 2-D numeric matrices \[LongDash] no \
toolbox imports.\n",
   "\[Bullet] ",
   StyleBox["DynamicModule", "Input"], " + ",
   StyleBox["PopupMenu", "Input"], " + ",
   StyleBox["Button", "Input"],
   " makes a cross-platform live downloader in thirty lines that \
works anywhere Wolfram is installed.\n",
   "\[Bullet] An ",
   StyleBox["asinh", "Input"],
   " stretch is the textbook recipe for turning raw telescope counts \
into a viewable RGB composite; three ",
   StyleBox["Reverse[..]", "Input"],
   " channels and a ",
   StyleBox["Transpose[_, {3, 1, 2}]", "Input"],
   " produce the final image in one line.\n",
   "\[Bullet] Images shipped with the notebook come from NASA \
SkyView; acknowledgements in \[Section]8."}], "Text"],

 Cell["10. References and further reading", "Section"],

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
the asinh-stretch colour-composite recipe used in \[Section]3.\n",
   "\[Bullet] This repository: ",
   ButtonBox["github.com/mthiel74/AstronomyImages",
     BaseStyle -> "Hyperlink",
     ButtonData -> {URL["https://github.com/mthiel74/AstronomyImages"],
                    None}],
   " \[LongDash] MIT-licensed code, bundled SkyView imagery."}],
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
