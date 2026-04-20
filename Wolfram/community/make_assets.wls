#!/usr/bin/env wolframscript
(* Regenerate the community-post assets (FITS, PNGs, stats JSON) using
   only Wolfram Language tooling. Mirrors what images/make_assets.py did
   on the Python side, but routed through HorseheadSkyView.wl so the
   notebook's code paths are exercised end-to-end.

   Usage:   wolframscript -file Wolfram/community/make_assets.wls
   Output:  Wolfram/community/images/{*.fits, *.png, analysis_stats.json}
*)

hereDir = DirectoryName[If[$ScriptCommandLine =!= {}, $InputFileName,
                           NotebookFileName[]]];
pkgDir  = FileNameJoin[{hereDir, ".."}];
imgDir  = FileNameJoin[{hereDir, "images"}];
If[!DirectoryQ[imgDir], CreateDirectory[imgDir]];

(* Load the Wolfram port but suppress its script-mode entry point by
   setting $Notebooks temporarily. *)
$Notebooks = True;
Get[FileNameJoin[{pkgDir, "HorseheadSkyView.wl"}]];

(* Robust FITS loader: cope with Association / nested-list returns. *)
loadFITSData[file_] := With[{raw = Quiet @ Import[file, {"FITS", "Data"}]},
  Which[
    AssociationQ[raw],                           First[Values[raw]],
    ListQ[raw] && Depth[raw] >= 4,               First[raw],
    True,                                        raw
  ]];

(* Render a PNG preview from a numeric image matrix. *)
savePreview[data_, path_, title_, cmap_: "GrayTones"] := Module[
  {flat, lo, hi},
  flat = DeleteCases[Flatten[data], _?(Not@*NumericQ)];
  {lo, hi} = Quantile[flat, {0.01, 0.995}];
  Export[path,
    ArrayPlot[Reverse[data],
      ColorFunction -> ColorData[cmap],
      PlotRange -> {lo, hi},
      Frame -> True,
      FrameLabel -> {"X (pixels)", "Y (pixels)"},
      PlotLabel -> Style[title, Bold, 12],
      ImageSize -> 640]];
  Print["  wrote ", path];
];

(* ---------------- 1. Real SkyView cutouts ---------------- *)

raHH  = 85.2479;
decHH = -2.4583;

fetchSurvey[survey_] := Module[
  {bytes, safe, fitsPath, data},
  safe = StringReplace[survey, {" " -> "_", "-" -> "_"}];
  fitsPath = FileNameJoin[{imgDir, "horsehead_" <> safe <> ".fits"}];
  Print["- trying SkyView survey ", survey];
  bytes = Quiet @ Catch[
    HorseheadSkyView`querySkyView[raHH, decHH,
      "Survey" -> survey, "Pixels" -> 600, "FOV" -> 0.5,
      "Timeout" -> 90],
    HorseheadSkyView`Private`querySkyView];
  If[!ByteArrayQ[bytes] || Length[bytes] < 8000,
    Print["  ! download failed for ", survey]; Return[$Failed]];
  BinaryWrite[fitsPath, bytes]; Close[fitsPath];
  data = loadFITSData[fitsPath];
  If[!MatrixQ[data, NumericQ],
    Print["  ! could not decode FITS: ", fitsPath]; Return[$Failed]];
  savePreview[data,
    FileNameJoin[{imgDir, "horsehead_" <> safe <> ".png"}],
    "Horsehead Nebula - " <> survey];
  data
];

realData = None;
Do[
  With[{d = fetchSurvey[s]},
    If[d =!= $Failed && realData === None, realData = d]],
  {s, {"DSS2 Red", "2MASS-J"}}];

(* ---------------- 2. Simulated Horsehead ---------------- *)

SeedRandom[42];
size = 600;
xs   = Subdivide[-1., 1., size - 1];
ys   = xs;
sim = Module[{X, Y, bg, mask, tau, img},
  X = Table[x, {y, ys}, {x, xs}];
  Y = Table[y, {y, ys}, {x, xs}];
  bg = 100. + 150. Exp[-(X^2 + Y^2)/0.5];
  Do[bg = bg + RandomReal[{10, 30}]
              Sin[RandomReal[{5, 15}]
                  (X Cos[RandomReal[{0, Pi}]] + Y Sin[RandomReal[{0, Pi}]])],
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

(* ---------------- 3. Analysis panel ---------------- *)

analysisSrc = If[realData =!= None, realData, sim];
flat = DeleteCases[Flatten[analysisSrc], _?(Not@*NumericQ)];
{lo, hi} = Quantile[flat, {0.01, 0.995}];
threshold = Quantile[flat, 0.30];
darkMask  = UnitStep[threshold - analysisSrc];

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
    Row[{"Threshold = ", NumberForm[N@threshold, {5, 1}]}]]}}, Spacings -> 1];

Export[FileNameJoin[{imgDir, "horsehead_analysis_panel.png"}], panel,
       ImageResolution -> 140];
Print["  wrote horsehead_analysis_panel.png"];

(* ---------------- 4. Stats JSON ---------------- *)

bgMean  = Mean[Pick[flat, Flatten[darkMask], 0]];
nebMean = Mean[Pick[flat, Flatten[darkMask], 1]];
tauEst  = -Log[(nebMean + 10.^-6)/(bgMean + 10.^-6)];
stats = <|
  "source"                    -> If[realData =!= None, "real", "simulated"],
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

Print["done; real data: ", realData =!= None];
