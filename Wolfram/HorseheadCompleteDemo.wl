(* ::Package:: *)

(* ::Title:: *)
(*HorseheadCompleteDemo.wl*)

(* ::Text:: *)
(*Wolfram Language translation of horsehead_complete_demo.py.*)
(*Simulates an H-alpha observation of the Horsehead Nebula, writes a FITS file*)
(*with a minimal WCS header, analyses the dark-nebula statistics, produces a*)
(*9-panel visualisation, and dumps API-call examples as JSON.*)


BeginPackage["HorseheadCompleteDemo`"];

createSimulatedHorseheadImage::usage =
  "createSimulatedHorseheadImage[size, noiseLevel] returns an Association with keys \
\"Image\", \"Background\" and \"OpticalDepth\" giving a physics-based simulation of the \
Horsehead dark nebula.";

createFITSWithWCS::usage =
  "createFITSWithWCS[data, ra, dec, pixelScale] returns an Association with FITS primary \
header cards and the image data, suitable for Export[file, ..., \"FITS\"].";

analyzeDarkNebula::usage =
  "analyzeDarkNebula[image, thresholdPercentile] prints dark-nebula diagnostics and \
returns {stats, darkMask}.";

plotComprehensiveAnalysis::usage =
  "plotComprehensiveAnalysis[image, background, opticalDepth, darkMask] returns a \
comprehensive 9-panel Graphics object.";

demonstrateAPICalls::usage =
  "demonstrateAPICalls[] prints sample SkyView and LCO request structures and returns \
them as an Association.";

horseheadCompleteDemoMain::usage =
  "horseheadCompleteDemoMain[outDir] reproduces the full Python demonstration: FITS + \
PNG + JSON outputs.";

Begin["`Private`"];

HORSEHEAD$RA  = 85.2479;
HORSEHEAD$DEC = -2.4583;
HORSEHEAD$SIZE$ARCMIN = 8;

(* ---------------- Simulation ---------------- *)

createSimulatedHorseheadImage[size_Integer:1000, noiseLevel_:0.02] := Module[
  {x, y, xy, X, Y, background, angles, freqs, amps, mask, tau, transmitted,
   image, readout, finite},

  Print["Simulating Horsehead Nebula observation..."];
  Print["  Image size: ", size, "x", size, " pixels"];
  Print["  Noise level: ", noiseLevel];

  x = Subdivide[-1., 1., size - 1];
  y = Subdivide[-1., 1., size - 1];
  X = ConstantArray[x, size];
  Y = Transpose[ConstantArray[y, size]];

  (* Background emission nebula (IC 434) - H-alpha emission *)
  background = 100. + 150. * Exp[-(X^2 + Y^2)/0.5];

  (* Filamentary structure *)
  angles = RandomReal[{0, Pi}, 5];
  freqs  = RandomReal[{5, 15}, 5];
  amps   = RandomReal[{10, 40}, 5];
  Do[
    background = background +
      amps[[i]] * Sin[freqs[[i]] * (X * Cos[angles[[i]]] + Y * Sin[angles[[i]]])],
    {i, 5}
  ];

  (* Horsehead mask - simplified *)
  mask = ConstantArray[0., {size, size}];
  (* Head *)
  mask = mask +
    Exp[-((X - 0.1)^2/0.3^2 + (Y + 0.2)^2/0.4^2) * 8];
  (* Neck *)
  mask = mask +
    0.7 * Exp[-((X - 0.15)^2/0.15^2 + (Y - 0.15)^2/0.35^2) * 6];

  tau         = 3. * mask;
  transmitted = Exp[-tau];
  image = background * transmitted;

  (* Poisson photon noise: approximate by Normal(mean, sqrt(mean)) to stay fast *)
  image = MapThread[Max[#1, 0] &, {
    image + Sqrt[Abs[image]] * RandomVariate[NormalDistribution[0., 1.], {size, size}],
    0. * image
  }, 2];

  readout = RandomVariate[NormalDistribution[0., noiseLevel * Mean[Flatten[image]]],
    {size, size}];
  image = image + readout;
  image = Map[Max[#, 0.] &, image, {2}];

  finite = Flatten[image];
  Print["  Intensity range: [", NumberForm[Min[finite], {5, 1}], ", ",
                                  NumberForm[Max[finite], {5, 1}], "] ADU"];
  Print["  Mean background: ", NumberForm[Mean[Flatten[background]], {5, 1}], " ADU"];
  Print["  Contrast ratio: ", NumberForm[Max[finite]/(Min[finite] + 10.^-6), {5, 2}]];

  <| "Image" -> image, "Background" -> background, "OpticalDepth" -> tau |>
];

(* ---------------- FITS WCS ---------------- *)

createFITSWithWCS[imageData_?MatrixQ, ra_?NumericQ, dec_?NumericQ,
                  pixelScale_:1.0] := Module[
  {ny, nx, header},
  {ny, nx} = Dimensions[imageData];
  header = <|
    "OBJECT"   -> "Horsehead Nebula",
    "TELESCOP" -> "Simulated Observatory",
    "INSTRUME" -> "Virtual CCD",
    "FILTER"   -> "H-alpha",
    "EXPTIME"  -> 300.0,
    "DATE-OBS" -> "2025-11-17T00:00:00",
    "CTYPE1"   -> "RA---TAN",
    "CTYPE2"   -> "DEC--TAN",
    "CRVAL1"   -> N[ra],
    "CRVAL2"   -> N[dec],
    "CRPIX1"   -> nx/2.,
    "CRPIX2"   -> ny/2.,
    "CDELT1"   -> -pixelScale/3600.,
    "CDELT2"   -> pixelScale/3600.,
    "CUNIT1"   -> "deg",
    "CUNIT2"   -> "deg",
    "EQUINOX"  -> 2000.0
  |>;
  <| "Data" -> imageData, "Header" -> header |>
];

exportFITS[fits_Association, filename_String] := Module[{tmp},
  (* Wolfram can Export a matrix to FITS; header cards are best-effort. *)
  Export[filename, fits["Data"], "FITS"];
  Print["\[Checkmark] FITS saved: ", filename];
  filename
];

(* ---------------- Analysis ---------------- *)

analyzeDarkNebula[image_?MatrixQ, thresholdPercentile_:30] := Module[
  {flat, threshold, darkMask, total, darkN, darkFrac, backMean, nebMean,
   ratio, tauEst, stats},
  Print["\n--- Dark Nebula Analysis ---"];

  flat      = Flatten[image];
  threshold = Quantile[flat, thresholdPercentile/100.];
  darkMask  = Map[# < threshold &, image, {2}];

  darkN    = Count[Flatten[darkMask], True];
  total    = Length[flat];
  darkFrac = darkN / total;

  backMean = Mean @ Pick[flat, Flatten[darkMask], False];
  nebMean  = Mean @ Pick[flat, Flatten[darkMask], True];

  ratio  = (nebMean + 10.^-6) / (backMean + 10.^-6);
  tauEst = -Log[ratio];

  stats = <|
    "dark_region_fraction" -> N[darkFrac],
    "background_intensity" -> N[backMean],
    "nebula_intensity"     -> N[nebMean],
    "contrast"             -> N[backMean/nebMean],
    "estimated_optical_depth" -> N[tauEst],
    "extinction_magnitudes" -> N[tauEst * 1.086]
  |>;

  Print["Dark region coverage: ", NumberForm[100 N[darkFrac], {5, 1}], "% of image"];
  Print["Background intensity: ", NumberForm[N[backMean], {5, 1}], " ADU"];
  Print["Nebula intensity:     ", NumberForm[N[nebMean], {5, 1}], " ADU"];
  Print["Intensity contrast:   ", NumberForm[N[stats["contrast"]], {5, 2}], "x"];
  Print["Estimated optical depth (tau): ", NumberForm[N[tauEst], {5, 2}]];
  Print["Estimated extinction: ", NumberForm[N[stats["extinction_magnitudes"]], {5, 2}],
        " magnitudes"];

  {stats, darkMask}
];

(* ---------------- Plotting ---------------- *)

arrayPlot[data_, opts___] := ArrayPlot[Reverse[data],
  ColorFunction -> "GrayTones",
  DataReversed  -> False,
  Frame         -> True,
  ImageSize     -> 320,
  opts];

plotComprehensiveAnalysis[image_, background_, tau_, darkMask_] := Module[
  {midRow, midCol, step, dims, panels, logImg, histPlot, xHist, hist, maskedView,
   surface3D, X3, Y3, Z3},

  dims = Dimensions[image];
  midRow = Floor[dims[[1]] / 2];
  midCol = Floor[dims[[2]] / 2];
  step = 20;

  logImg = Log10[image + 1];

  hist = Histogram[Flatten[image], 100,
    ChartStyle -> RGBColor[0.27, 0.51, 0.71],
    GridLines  -> {{{Quantile[Flatten[image], 0.30], Red}}, None},
    Frame      -> True,
    FrameLabel -> {"Pixel Intensity (ADU)", "Frequency"},
    PlotLabel  -> "Intensity Distribution",
    ImageSize  -> 320
  ];

  maskedView = MapThread[If[#2, #1, Indeterminate] &, {image, darkMask}, 2];

  (* 3D surface: downsample *)
  Z3 = image[[;; ;; step, ;; ;; step]];
  X3 = Range[0, dims[[2]] - 1, step];
  Y3 = Range[0, dims[[1]] - 1, step];
  surface3D = ListPlot3D[Z3,
    ColorFunction -> "TerrainColors",
    Mesh -> None,
    ImageSize -> 320,
    AxesLabel -> {"X", "Y", "Intensity"},
    PlotLabel -> "3D Intensity Surface",
    ViewPoint -> {1.5, -2.0, 1.2}
  ];

  panels = {
    { Labeled[arrayPlot[image, ColorFunction -> "GrayTones"],
        "Simulated Observation\n(H-alpha filter)", Top],
      Labeled[arrayPlot[background, ColorFunction -> "SunsetColors"],
        "Background Emission Nebula\n(IC 434, intrinsic)", Top],
      Labeled[arrayPlot[tau, ColorFunction -> "SouthwestColors"],
        "Dust Optical Depth\n(\[Tau], absorption)", Top]
    },
    { Labeled[arrayPlot[logImg, ColorFunction -> "GrayTones"],
        "Log-Scaled Observation\n(enhanced contrast)", Top],
      Labeled[arrayPlot[maskedView, ColorFunction -> "RedBlueTones"],
        "Segmented Dark Nebula\n(absorption regions)", Top],
      hist
    },
    { ListLinePlot[{image[[midRow]], background[[midRow]]},
        PlotLegends -> {"Observed", "Background"},
        Frame -> True, GridLines -> Automatic,
        FrameLabel -> {"Pixel X", "Intensity (ADU)"},
        PlotLabel -> "Horizontal Cross-Section (Y=" <> ToString[midRow] <> ")",
        ImageSize -> 320],
      ListLinePlot[{image[[All, midCol]], background[[All, midCol]]},
        PlotLegends -> {"Observed", "Background"},
        Frame -> True, GridLines -> Automatic,
        FrameLabel -> {"Pixel Y", "Intensity (ADU)"},
        PlotLabel -> "Vertical Cross-Section (X=" <> ToString[midCol] <> ")",
        ImageSize -> 320],
      surface3D
    }
  };

  Labeled[
    Grid[panels, Spacings -> {1, 1}],
    Style["Horsehead Nebula (Barnard 33) - Complete Analysis\nSimulated H-alpha Observation",
      Bold, 14],
    Top
  ]
];

(* ---------------- API call examples ---------------- *)

demonstrateAPICalls[] := Module[{examples},
  Print["\n", StringJoin[ConstantArray["=", 70]]];
  Print["EXAMPLE API CALLS FOR REAL SERVICES"];
  Print[StringJoin[ConstantArray["=", 70]]];

  examples = <|
    "SkyView" -> <|
      "url" -> "https://skyview.gsfc.nasa.gov/cgi-bin/images",
      "method" -> "GET",
      "params" -> <|
        "Position" -> ToString[HORSEHEAD$RA] <> "," <> ToString[HORSEHEAD$DEC],
        "Survey"   -> "DSS2 Red",
        "Pixels"   -> "1000,1000",
        "Size"     -> 0.5,
        "Return"   -> "FITS"
      |>,
      "wolfram" -> "
fits = Import[
  \"https://skyview.gsfc.nasa.gov/cgi-bin/images?\" <>
  URLQueryEncode[<|
    \"Position\" -> \"85.2479,-2.4583\",
    \"Survey\"   -> \"DSS2 Red\",
    \"Pixels\"   -> \"1000,1000\",
    \"Size\"     -> 0.5,
    \"Return\"   -> \"FITS\"
  |>],
  \"FITS\"];
"
    |>,
    "Las Cumbres Observatory" -> <|
      "url" -> "https://observe.lco.global/api/requestgroups/",
      "method" -> "POST",
      "authentication" -> "Token YOUR_API_TOKEN",
      "body" -> <|
        "name" -> "Horsehead Nebula Observation",
        "proposal" -> "EDUCATION_PROJECT",
        "ipp_value" -> 1.0,
        "operator" -> "SINGLE",
        "observation_type" -> "NORMAL",
        "requests" -> {<|
          "target" -> <|
            "name" -> "Horsehead Nebula",
            "type" -> "ICRS",
            "ra"   -> HORSEHEAD$RA,
            "dec"  -> HORSEHEAD$DEC
          |>,
          "configurations" -> {<|
            "type" -> "EXPOSE",
            "instrument_type" -> "1M0-SCICAM-SINISTRO",
            "exposure_time" -> 300,
            "exposure_count" -> 3,
            "filter" -> "rp"
          |>},
          "windows" -> {<|
            "start" -> "2025-11-17T00:00:00Z",
            "end"   -> "2025-11-18T00:00:00Z"
          |>},
          "location" -> <|"telescope_class" -> "1m0"|>
        |>}
      |>,
      "wolfram" -> "
response = URLRead[HTTPRequest[
  \"https://observe.lco.global/api/requestgroups/\",
  <|
    \"Method\"      -> \"POST\",
    \"Headers\"     -> {\"Authorization\" -> \"Token YOUR_API_TOKEN\"},
    \"ContentType\" -> \"application/json\",
    \"Body\"        -> ExportString[observationRequest, \"JSON\"]
  |>]];
"
    |>
  |>;

  KeyValueMap[
    Function[{service, details},
      Print["\n--- ", service, " ---"];
      Print["Endpoint: ", details["url"]];
      Print["Method:   ", details["method"]];
      If[ KeyExistsQ[details, "authentication"],
        Print["Auth:     ", details["authentication"]]
      ];
      Print["\nWolfram Example:"];
      Print[details["wolfram"]];
    ],
    examples
  ];
  examples
];

(* ---------------- Main ---------------- *)

horseheadCompleteDemoMain[outDir_String:"."] := Module[
  {sim, image, background, tau, fitsObj, fitsPath, stats, darkMask, statsPath,
   fig, plotPath, apiExamples, apiDocPath},

  Print[StringJoin[ConstantArray["=", 70]]];
  Print["REMOTE TELESCOPE & VIRTUAL OBSERVATORY DEMONSTRATION"];
  Print[StringJoin[ConstantArray["=", 70]]];
  Print["\nTarget: Horsehead Nebula (Barnard 33)"];
  Print["Coordinates: RA=", HORSEHEAD$RA, "\[Degree] (",
        NumberForm[HORSEHEAD$RA/15., {3, 2}], "h), Dec=", HORSEHEAD$DEC, "\[Degree]"];
  Print["Constellation: Orion"];
  Print["Type: Dark nebula (dust absorption feature)"];
  Print["Size: ~", HORSEHEAD$SIZE$ARCMIN, "' x 5'"];

  sim        = createSimulatedHorseheadImage[800];
  image      = sim["Image"];
  background = sim["Background"];
  tau        = sim["OpticalDepth"];

  Print["\n--- Creating FITS File ---"];
  fitsObj  = createFITSWithWCS[image, HORSEHEAD$RA, HORSEHEAD$DEC, 1.0];
  fitsPath = FileNameJoin[{outDir, "horsehead_simulated.fits"}];
  exportFITS[fitsObj, fitsPath];

  {stats, darkMask} = analyzeDarkNebula[image];
  statsPath = FileNameJoin[{outDir, "horsehead_analysis.json"}];
  Export[statsPath, stats, "JSON"];
  Print["\[Checkmark] Statistics saved: ", statsPath];

  Print["\n--- Creating Visualizations ---"];
  fig = plotComprehensiveAnalysis[image, background, tau, darkMask];
  plotPath = FileNameJoin[{outDir, "horsehead_analysis.png"}];
  Export[plotPath, Rasterize[fig, ImageResolution -> 150]];
  Print["\[Checkmark] Analysis plot saved: ", plotPath];

  apiExamples = demonstrateAPICalls[];
  apiDocPath  = FileNameJoin[{outDir, "api_examples.json"}];
  Export[apiDocPath, apiExamples, "JSON"];
  Print["\n\[Checkmark] API documentation saved: ", apiDocPath];

  Print["\n", StringJoin[ConstantArray["=", 70]]];
  Print["DEMONSTRATION COMPLETE"];
  Print[StringJoin[ConstantArray["=", 70]]];
  Print["\nFiles created:"];
  Print["  1. ", fitsPath, " - FITS image"];
  Print["  2. ", plotPath, " - Comprehensive analysis"];
  Print["  3. ", statsPath, " - Quantitative measurements"];
  Print["  4. ", apiDocPath, " - API call examples"];

  Print["\nNext steps to use REAL telescopes:"];
  Print["  1. SkyView (free): see HorseheadSkyView.wl"];
  Print["  2. LCO (free for education): Sign up at observe.lco.global"];
  Print["  3. iTelescope (paid): Professional-grade remote access"];
];

End[];
EndPackage[];

If[ !$Notebooks, HorseheadCompleteDemo`horseheadCompleteDemoMain[DirectoryName[$InputFileName]] ]
