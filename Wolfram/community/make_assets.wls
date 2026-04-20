#!/usr/bin/env wolframscript
(* Regenerate the community-post assets using only Wolfram Language.
   This version builds:

     1. Horsehead Nebula in DSS2 Red + 2MASS-J + simulated        (core)
     2. Multi-wavelength Horsehead panel: DSS2 Red, 2MASS-J,
        2MASS-K, WISE 3.4 \[LongDash] one object, four bands
     3. 2MASS J/H/K Horsehead RGB composite
     4. Messier gallery: M42, M45, M31, M16, M51, M1 in DSS2 Red
     5. Analysis panel + stats JSON

   Downloads are idempotent: if the FITS already exists we skip the
   network round-trip.

   Usage:   wolframscript -file Wolfram/community/make_assets.wls
*)

hereDir = DirectoryName[If[$ScriptCommandLine =!= {}, $InputFileName,
                           NotebookFileName[]]];
pkgDir  = FileNameJoin[{hereDir, ".."}];
imgDir  = FileNameJoin[{hereDir, "images"}];
If[!DirectoryQ[imgDir], CreateDirectory[imgDir]];

$Notebooks = True;   (* suppress script-mode main in the package *)
Get[FileNameJoin[{pkgDir, "HorseheadSkyView.wl"}]];

loadFITSData[file_] := With[{raw = Quiet @ Import[file, {"FITS", "Data"}]},
  Which[
    AssociationQ[raw],                           First[Values[raw]],
    ListQ[raw] && Depth[raw] >= 4,               First[raw],
    True,                                        raw]];

safeName[s_] := StringReplace[s, {" " -> "_", "-" -> "_", "." -> "p",
                                   "(" -> "", ")" -> ""}];

(* ---------- core SkyView fetcher with filesystem caching ---------- *)

fetch[ra_, dec_, survey_, fov_, pixels_, tag_String] := Module[
  {safe, fitsPath, bytes, data},
  safe     = safeName[tag] <> "_" <> safeName[survey];
  fitsPath = FileNameJoin[{imgDir, safe <> ".fits"}];
  If[FileExistsQ[fitsPath] && FileByteCount[fitsPath] > 8000,
    Print["  cache hit: ", FileNameTake[fitsPath]];
    data = loadFITSData[fitsPath];
    Return[{fitsPath, data}]];
  Print["  downloading ", tag, " / ", survey, " (",
        fov, "\[Degree], ", pixels, "px)"];
  bytes = Quiet @ Catch[
    HorseheadSkyView`querySkyView[ra, dec,
      "Survey" -> survey, "Pixels" -> pixels, "FOV" -> fov,
      "Timeout" -> 120],
    HorseheadSkyView`Private`querySkyView];
  If[!ByteArrayQ[bytes] || Length[bytes] < 8000,
    Print["  ! failed: ", tag, " / ", survey]; Return[{$Failed, $Failed}]];
  BinaryWrite[fitsPath, bytes]; Close[fitsPath];
  data = loadFITSData[fitsPath];
  If[!MatrixQ[data, NumericQ],
    Print["  ! could not decode FITS: ", fitsPath]; Return[{fitsPath, $Failed}]];
  {fitsPath, data}
];

(* ---------- rendering helpers ---------- *)

percClip[data_, lo_: 0.01, hi_: 0.995] := Module[{flat},
  flat = DeleteCases[Flatten[data], _?(Not@*NumericQ)];
  Quantile[flat, {lo, hi}]];

savePreview[data_, path_, title_, cmap_: "GrayTones"] := Module[
  {lohi},
  lohi = percClip[data];
  Export[path,
    ArrayPlot[Reverse[data],
      ColorFunction -> ColorData[cmap],
      PlotRange -> lohi,
      Frame -> True,
      FrameLabel -> {"X (pixels)", "Y (pixels)"},
      PlotLabel -> Style[title, Bold, 12],
      ImageSize -> 640]];
  Print["  wrote ", FileNameTake[path]]];

(* Normalise a 2-D array to [0,1] using asinh stretch on percentile-
   clipped values \[LongDash] the textbook trick for turning raw
   telescope counts into visually balanced RGB channels. *)
asinhStretch[data_, quantileLo_: 0.5, quantileHi_: 0.998,
             softening_: 0.02] := Module[
  {flat, lo, hi, norm, stretched, sMin, sMax},
  flat = DeleteCases[Flatten[data], _?(Not@*NumericQ)];
  {lo, hi} = Quantile[flat, {quantileLo, quantileHi}];
  norm = Clip[(data - lo)/(hi - lo + 10.^-12), {0, 1}];
  stretched = ArcSinh[norm / softening];
  {sMin, sMax} = MinMax[Flatten[stretched]];
  Clip[(stretched - sMin)/(sMax - sMin + 10.^-12), {0, 1}]];

(* =================================================================
   1. Horsehead Nebula: core panels
   ================================================================= *)

raHH  = 85.2479;
decHH = -2.4583;

Print["---- Core Horsehead panels"];
{ , dataDSSRed}  = fetch[raHH, decHH, "DSS2 Red",  0.5, 600, "horsehead"];
{ , dataMASSJ}   = fetch[raHH, decHH, "2MASS-J",   0.5, 600, "horsehead"];

If[MatrixQ[dataDSSRed, NumericQ],
  savePreview[dataDSSRed,
    FileNameJoin[{imgDir, "horsehead_DSS2_Red.png"}],
    "Horsehead Nebula - DSS2 Red"]];
If[MatrixQ[dataMASSJ, NumericQ],
  savePreview[dataMASSJ,
    FileNameJoin[{imgDir, "horsehead_2MASS_J.png"}],
    "Horsehead Nebula - 2MASS J (1.25\[Mu]m)"]];

(* =================================================================
   2. Multi-wavelength Horsehead (DSS Red, 2MASS J, K, WISE 3.4)
   ================================================================= *)

Print["---- Multi-wavelength Horsehead"];
{ , dataMASSK}  = fetch[raHH, decHH, "2MASS-K",  0.5, 600, "horsehead"];
{ , dataWISE34} = fetch[raHH, decHH, "WISE 3.4", 0.5, 600, "horsehead"];

multiWavePanel = Module[{panels, mk},
  mk[data_, label_, cmap_] :=
    If[MatrixQ[data, NumericQ],
      With[{lohi = percClip[data]},
        Labeled[
          ArrayPlot[Reverse[data],
            ColorFunction -> ColorData[cmap],
            PlotRange -> lohi, Frame -> False, ImageSize -> 380],
          Style[label, Bold, 13]]],
      Framed[Style["[missing]", Red], ImageSize -> {380, 380}]];
  panels = {
    mk[dataDSSRed,  "DSS2 Red (656 nm)",          "GrayTones"],
    mk[dataMASSJ,   "2MASS J (1.25 \[Mu]m)",      "GrayTones"],
    mk[dataMASSK,   "2MASS K (2.17 \[Mu]m)",      "SunsetColors"],
    mk[dataWISE34,  "WISE W1 (3.4 \[Mu]m)",       "SunsetColors"]};
  Labeled[
    Grid[Partition[panels, 2], Spacings -> 1],
    Style["Horsehead Nebula: one object, four wavelengths",
          Bold, 14], Top]];

Export[FileNameJoin[{imgDir, "horsehead_multiwavelength.png"}],
  multiWavePanel, ImageResolution -> 140];
Print["  wrote horsehead_multiwavelength.png"];

(* =================================================================
   3. 2MASS J/H/K RGB composite (using asinh stretch)
   ================================================================= *)

Print["---- 2MASS J/H/K RGB composite"];
{ , dataMASSH} = fetch[raHH, decHH, "2MASS-H", 0.5, 600, "horsehead"];

If[AllTrue[{dataMASSJ, dataMASSH, dataMASSK}, MatrixQ[#, NumericQ] &],
  Module[{r, g, b, rgb},
    r = asinhStretch[dataMASSK];    (* K \[Rule] red   *)
    g = asinhStretch[dataMASSH];    (* H \[Rule] green *)
    b = asinhStretch[dataMASSJ];    (* J \[Rule] blue  *)
    rgb = Transpose[Map[Reverse, {r, g, b}], {3, 1, 2}];
    Export[FileNameJoin[{imgDir, "horsehead_RGB_JHK.png"}],
      Rasterize[
        Labeled[
          Image[rgb, Interleaving -> True, ImageSize -> 640],
          Style["Horsehead Nebula: 2MASS J/H/K false-colour composite",
                Bold, 13], Top],
        ImageSize -> 720]];
    Print["  wrote horsehead_RGB_JHK.png"]],
  Print["  ! RGB skipped, missing 2MASS band"]];

(* =================================================================
   4. Messier gallery
   ================================================================= *)

Print["---- Messier gallery"];

gallery = <|
  "M42 (Orion Nebula)"    -> <|"ra" -> 83.8213, "dec" -> -5.3903,
                              "fov" -> 1.0,  "cmap" -> "Rainbow",
                              "caption" -> "Stellar nursery, 1344 ly away"|>,
  "M45 (Pleiades)"        -> <|"ra" -> 56.870,  "dec" -> 24.105,
                              "fov" -> 2.5,  "cmap" -> "DeepSeaColors",
                              "caption" -> "Open cluster, ~110 Myr old"|>,
  "M31 (Andromeda)"       -> <|"ra" -> 10.685,  "dec" -> 41.269,
                              "fov" -> 3.0,  "cmap" -> "SunsetColors",
                              "caption" -> "Nearest spiral galaxy, 2.5 Mly"|>,
  "M16 (Eagle Nebula)"    -> <|"ra" -> 274.700, "dec" -> -13.807,
                              "fov" -> 0.5,  "cmap" -> "TemperatureMap",
                              "caption" -> "Pillars of Creation, HST target"|>,
  "M51 (Whirlpool)"       -> <|"ra" -> 202.469, "dec" -> 47.196,
                              "fov" -> 0.2,  "cmap" -> "AuroraColors",
                              "caption" -> "Grand-design spiral, 23 Mly"|>,
  "M1 (Crab Nebula)"      -> <|"ra" -> 83.633,  "dec" -> 22.014,
                              "fov" -> 0.15, "cmap" -> "SolarColors",
                              "caption" -> "Supernova remnant of 1054 AD"|>
|>;

galleryPanels = {};

Do[
  Module[{name = k, meta = gallery[k], d, fpath, lohi, panel},
    Print["- ", name];
    {fpath, d} = fetch[meta["ra"], meta["dec"], "DSS2 Red",
                       meta["fov"], 600, name];
    If[MatrixQ[d, NumericQ],
      lohi = percClip[d];
      panel = Labeled[
        ArrayPlot[Reverse[d],
          ColorFunction -> ColorData[meta["cmap"]],
          PlotRange -> lohi,
          Frame -> False, ImageSize -> 360],
        Column[{
          Style[name, Bold, 12],
          Style[meta["caption"], Italic, FontSize -> 10,
                FontColor -> Gray]},
          Alignment -> Center]];
      AppendTo[galleryPanels, panel];
      (* Also write an individual PNG for the one-off gallery image *)
      Export[FileNameJoin[{imgDir,
          "gallery_" <> safeName[name] <> ".png"}],
        ArrayPlot[Reverse[d],
          ColorFunction -> ColorData[meta["cmap"]],
          PlotRange -> lohi,
          Frame -> False, ImageSize -> 640]];
      Print["  wrote individual gallery_", safeName[name], ".png"]]],
  {k, Keys[gallery]}];

If[Length[galleryPanels] >= 4,
  Export[FileNameJoin[{imgDir, "messier_gallery.png"}],
    Labeled[
      Grid[Partition[galleryPanels, 3, 3, 1, ""],
           Spacings -> {1.2, 1.5}],
      Style["A six-target tour through NASA SkyView \
(DSS2 Red; credit: STScI, see notebook \[Section]10)", Bold, 14], Top],
    ImageResolution -> 140];
  Print["  wrote messier_gallery.png (",
        Length[galleryPanels], " targets)"]];

(* =================================================================
   5. Simulated Horsehead, analysis panel, stats JSON
   ================================================================= *)

Print["---- Simulated Horsehead + analysis panel"];

SeedRandom[42];
size = 600;
xs   = Subdivide[-1., 1., size - 1];
sim = Module[{X, Y, bg, mask, tau, img},
  X = Table[x, {y, xs}, {x, xs}];
  Y = Table[y, {y, xs}, {x, xs}];
  bg = 100. + 150. Exp[-(X^2 + Y^2)/0.5];
  Do[bg = bg + RandomReal[{10, 30}]
              Sin[RandomReal[{5, 15}]
                  (X Cos[RandomReal[{0, Pi}]]
                 + Y Sin[RandomReal[{0, Pi}]])],
     {5}];
  mask = Exp[-(((X - 0.10)/0.30)^2 + ((Y + 0.20)/0.40)^2) 8.]
       + 0.7 Exp[-(((X - 0.15)/0.15)^2 + ((Y - 0.15)/0.35)^2) 6.];
  tau = 3.0 mask;
  img = bg Exp[-tau];
  img = Map[RandomVariate[PoissonDistribution[Max[#, 0.]]] &,
            img, {2}];
  img + RandomVariate[NormalDistribution[0., 0.02 Mean[Flatten[img]]],
                      {size, size}]];
sim = Clip[sim, {0., Infinity}];
Export[FileNameJoin[{imgDir, "horsehead_simulated.mx"}], sim];
savePreview[sim,
  FileNameJoin[{imgDir, "horsehead_simulated.png"}],
  "Simulated Horsehead Nebula"];

analysisSrc = If[MatrixQ[dataDSSRed, NumericQ], dataDSSRed, sim];
flat = DeleteCases[Flatten[analysisSrc], _?(Not@*NumericQ)];
{lo, hi}    = Quantile[flat, {0.01, 0.995}];
threshold   = Quantile[flat, 0.30];
darkMask    = UnitStep[threshold - analysisSrc];

panel = Grid[{
  {Labeled[ArrayPlot[Reverse[analysisSrc],
      ColorFunction -> ColorData["GrayTones"],
      PlotRange -> {lo, hi}, ImageSize -> 360, Frame -> False],
    "Raw image"],
   Labeled[ArrayPlot[Reverse[Log10[analysisSrc + 1]],
      ColorFunction -> "SunsetColors",
      ImageSize -> 360, Frame -> False],
    "log10 scaling"]},
  {Labeled[Show[
      ArrayPlot[Reverse[analysisSrc],
        ColorFunction -> ColorData["GrayTones"],
        PlotRange -> {lo, hi}, ImageSize -> 360, Frame -> False],
      ArrayPlot[Reverse[darkMask],
        ColorRules -> {0 -> Opacity[0],
                       1 -> RGBColor[1, 0, 0, 0.45]},
        ImageSize -> 360, Frame -> False]],
    "Dark-nebula segmentation (bottom 30%)"],
   Labeled[Histogram[flat, 80,
      ChartStyle -> RGBColor[0.27, 0.51, 0.71],
      GridLines -> {{{threshold, Directive[Red, Dashed]}}, None},
      Frame -> True,
      FrameLabel -> {"Pixel intensity (ADU)", "Count"},
      ImageSize -> 360],
    Row[{"Threshold = ", NumberForm[N@threshold, {5, 1}]}]]}},
  Spacings -> 1];

Export[FileNameJoin[{imgDir, "horsehead_analysis_panel.png"}],
  panel, ImageResolution -> 140];
Print["  wrote horsehead_analysis_panel.png"];

bgMean  = Mean[Pick[flat, Flatten[darkMask], 0]];
nebMean = Mean[Pick[flat, Flatten[darkMask], 1]];
tauEst  = -Log[(nebMean + 10.^-6)/(bgMean + 10.^-6)];
stats = <|
  "source"                    -> If[MatrixQ[dataDSSRed, NumericQ],
                                    "real", "simulated"],
  "shape"                     -> Dimensions[analysisSrc],
  "min"                       -> N@Min[flat],
  "max"                       -> N@Max[flat],
  "mean"                      -> N@Mean[flat],
  "median"                    -> N@Median[flat],
  "threshold_p30"             -> N@threshold,
  "background_mean"           -> N@bgMean,
  "dark_mean"                 -> N@nebMean,
  "estimated_optical_depth"   -> N@tauEst,
  "estimated_extinction_mag"  -> N[1.086 tauEst]
|>;
Export[FileNameJoin[{imgDir, "analysis_stats.json"}], stats, "JSON"];
Print["  wrote analysis_stats.json"];

(* =================================================================
   6. GUI preview \[LongDash] a static screenshot of the live
      downloader, for inclusion in the notebook's \[Section]8.2
      (since a Dynamic widget itself cannot be inlined as a figure).
   ================================================================= *)

Print["---- GUI preview"];

guiPreviewSurveys  = {"DSS2 Red", "DSS2 Blue", "DSS2 IR",
                      "2MASS-J", "2MASS-H", "2MASS-K",
                      "WISE 3.4", "WISE 4.6", "WISE 12", "WISE 22",
                      "GALEX Near UV"};
guiPreviewTargets  = {"Horsehead Nebula", "M42 Orion Nebula",
                      "M45 Pleiades", "M31 Andromeda",
                      "M16 Eagle Nebula", "M51 Whirlpool",
                      "M1 Crab Nebula", "M104 Sombrero",
                      "NGC 2024 Flame", "Veil Nebula", "(custom)"};
guiPreviewCmaps    = {"GrayTones", "SunsetColors", "TemperatureMap",
                      "Rainbow", "DeepSeaColors", "AuroraColors",
                      "SolarColors"};

(* Use the already-downloaded M42 cutout as the "result" shown on the
   right-hand plot panel. Falls back to the Horsehead DSS Red if M42
   isn't around. *)
guiPreviewFile =
  FirstCase[{FileNameJoin[{imgDir, "M42_Orion_Nebula_DSS2_Red.fits"}],
             FileNameJoin[{imgDir, "horsehead_DSS2_Red.fits"}]},
    f_String /; FileExistsQ[f] :> f, None];

If[guiPreviewFile =!= None,
  Module[{gd, gLo, gHi, controls, plot, gui},
    gd = loadFITSData[guiPreviewFile];
    {gLo, gHi} = percClip[gd];
    controls = Panel[Column[{
      Style["Live SkyView downloader", Bold, 14],
      Row[{"Target: ",  PopupMenu["M42 Orion Nebula",
          guiPreviewTargets, FieldSize -> 22]}],
      Row[{"RA (\[Degree]): ",  InputField[83.8213, Number, FieldSize -> 10]}],
      Row[{"Dec (\[Degree]): ", InputField[-5.3903, Number, FieldSize -> 10]}],
      Row[{"FOV (\[Degree]): ", Slider[1.0, {0.02, 5.0, 0.02}],
           " ", Style["1.00", FontSize -> 11]}],
      Row[{"Survey: ",      PopupMenu["DSS2 Red",  guiPreviewSurveys,
                                      FieldSize -> 22]}],
      Row[{"Pixels: ",      PopupMenu[600, {300, 500, 600, 800, 1000}]}],
      Row[{"Colour map: ",  PopupMenu["SunsetColors", guiPreviewCmaps]}],
      Row[{
        Button[Style["Download & plot", Bold, White], Null,
          Background -> RGBColor[0.12, 0.45, 0.80], ImageSize -> 200],
        Spacer[6],
        Button["Save FITS", Null, ImageSize -> 120]}],
      Style["\[Checkmark] {600, 600} image from DSS2 Red (653,760 bytes)",
            Italic, FontSize -> 11]
    }, Spacings -> 0.8], ImageSize -> 400];
    plot = ArrayPlot[Reverse[gd],
      ColorFunction -> ColorData["SunsetColors"],
      PlotRange -> {gLo, gHi},
      Frame -> True,
      PlotLabel -> Style["M42 Orion Nebula \[LongDash] DSS2 Red",
                         Bold, 12],
      FrameLabel -> {"X (pixels)", "Y (pixels)"},
      ImageSize -> 540];
    gui = Row[{controls, Spacer[12], plot}];
    Export[FileNameJoin[{imgDir, "gui_preview.png"}],
      Rasterize[gui, ImageResolution -> 140]];
    Print["  wrote gui_preview.png"]],
  Print["  ! no cutout available for GUI preview"]];

(* =================================================================
   7. Click-map preview \[LongDash] static shot of the \[Section]6.3
      click-anywhere sky map + the cutout that would result from
      a click near the Horsehead.
   ================================================================= *)

Print["---- Click-map preview"];

skyLandmarksPreview = <|
  "Horsehead"     -> {85.25,  -2.46},
  "M42 Orion"     -> {83.82,  -5.39},
  "M45 Pleiades"  -> {56.87,  24.11},
  "M31 Andromeda" -> {10.68,  41.27},
  "M16 Eagle"     -> {274.70, -13.81},
  "M51 Whirlpool" -> {202.47, 47.20},
  "M1 Crab"       -> {83.63,  22.01}|>;

skyMapPreview[selRA_, selDec_] := Graphics[{
    {Purple, Opacity[0.15], Thickness[0.03],
      Line[{#, -29 + 60 Sin[(# - 266) Degree]} & /@ Subdivide[0., 360., 99]]},
    {Blue, Opacity[0.4], Line[{{0, 0}, {360, 0}}]},
    {Orange, Opacity[0.4], Dashed,
      Line[{#, 23.4 Sin[# Degree]} & /@ Subdivide[0., 360., 99]]},
    KeyValueMap[
      Function[{n, p},
        {Red, PointSize[0.012], Point[p],
         Text[Style[n, 9, Darker[Red]], p, {0, -1.5}]}],
      skyLandmarksPreview],
    {Green, PointSize[0.02], Point[{selRA, selDec}],
     Thickness[0.004], Circle[{selRA, selDec}, 6]}},
  PlotRange -> {{360, 0}, {-90, 90}},
  Frame -> True,
  FrameLabel -> {"Right Ascension (\[Degree])", "Declination (\[Degree])"},
  AspectRatio -> 1/2,
  GridLines -> Automatic,
  GridLinesStyle -> Directive[Gray, Opacity[0.3], Dashed],
  ImageSize -> 680];

If[guiPreviewFile =!= None,
  Module[{gd, gLo, gHi, map, plot, status, preview},
    gd = loadFITSData[
      FileNameJoin[{imgDir, "horsehead_DSS2_Red.fits"}]];
    {gLo, gHi} = percClip[gd];
    map = skyMapPreview[85.25, -2.46];
    plot = ArrayPlot[Reverse[gd],
      ColorFunction -> ColorData["GrayTones"],
      PlotRange -> {gLo, gHi},
      Frame -> True,
      PlotLabel -> Style[
        "DSS2 Red @ RA=85.25\[Degree], Dec=\[Minus]2.46\[Degree]",
        Bold, 11],
      ImageSize -> 540];
    status = Style[
      "\[Checkmark] {600, 600} image from DSS2 Red", Italic, 11];
    preview = Column[{
      map,
      Row[{"RA: 85.250\[Degree]   Dec: \[Minus]2.460\[Degree]   \
FOV: 0.50\[Degree]   Survey: ", PopupMenu["DSS2 Red",
        {"DSS2 Red", "2MASS-J", "2MASS-K", "WISE 3.4"}], "   ",
        Button[Style["Download", Bold, White], Null,
          Background -> RGBColor[0.12, 0.45, 0.80], ImageSize -> 120]}],
      status,
      plot}, Spacings -> 1];
    Export[FileNameJoin[{imgDir, "click_map_preview.png"}],
      Rasterize[preview, ImageResolution -> 140]];
    Print["  wrote click_map_preview.png"]]];

Print["\nAll done."];
