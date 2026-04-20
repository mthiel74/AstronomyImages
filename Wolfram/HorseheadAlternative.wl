(* ::Package:: *)

(* ::Title:: *)
(*HorseheadAlternative.wl*)

(* ::Text:: *)
(*Wolfram Language translation of horsehead_alternative.py.*)
(*Queries the SDSS DR18 SkyServer ImgCutout service and the STScI DSS2 service for*)
(*JPEG renderings of the Horsehead Nebula region, then analyses the image intensity.*)


BeginPackage["HorseheadAlternative`"];

querySDSSImage::usage =
  "querySDSSImage[ra, dec, opts] returns a JPEG Image from SDSS SkyServer or $Failed. \
Options: \"Scale\" (arcsec/pix), \"Width\" (pix), \"Height\" (pix).";

queryDSSViaAlternate::usage =
  "queryDSSViaAlternate[ra, dec, sizeArcmin] returns a JPEG Image from STScI DSS2 \
or $Failed.";

displayAndAnalyze::usage =
  "displayAndAnalyze[img, title, savePath] saves side-by-side and histogram PNGs \
and prints basic intensity statistics.";

horseheadAlternativeMain::usage =
  "horseheadAlternativeMain[outDir] reproduces the Python demo using SDSS and DSS2.";

Begin["`Private`"];

HORSEHEAD$RA  = 85.2479;
HORSEHEAD$DEC = -2.4583;

Options[querySDSSImage] = {
  "Scale"  -> 0.4,
  "Width"  -> 512,
  "Height" -> 512,
  "Opt"    -> "G",
  "Timeout" -> 30
};

querySDSSImage[ra_?NumericQ, dec_?NumericQ, OptionsPattern[]] := Module[
  {baseURL, params, url, scale, width, height},
  scale  = OptionValue["Scale"];
  width  = OptionValue["Width"];
  height = OptionValue["Height"];
  baseURL = "http://skyserver.sdss.org/dr18/SkyServerWS/ImgCutout/getjpeg";
  params = <|
    "ra"     -> ToString[ra],
    "dec"    -> ToString[dec],
    "scale"  -> ToString[scale],
    "width"  -> ToString[width],
    "height" -> ToString[height],
    "opt"    -> OptionValue["Opt"]
  |>;
  url = baseURL <> "?" <> URLQueryEncode[params];

  Print["Querying SDSS DR18 for coordinates: RA=", ra, "\[Degree], Dec=", dec, "\[Degree]"];
  Print["Scale: ", scale, " arcsec/pixel, Size: ", width, "x", height, " pixels"];
  Print["Field of view: ~", NumberForm[width*scale/3600., {3, 2}], "\[Degree] x ",
        NumberForm[height*scale/3600., {3, 2}], "\[Degree]"];

  Module[{img},
    img = Quiet @ TimeConstrained[
      Import[url, "JPEG"], OptionValue["Timeout"], $Failed];
    If[ Head[img] === Image,
      Print["\[Checkmark] Successfully retrieved image"]; img,
      Print["\[Times] SDSS query failed"]; $Failed
    ]
  ]
];

queryDSSViaAlternate[ra_?NumericQ, dec_?NumericQ, sizeArcmin_:30] := Module[
  {baseURL, raHours, raH, raM, raS, decSign, decAbs, decD, decM, decS, params, url, img},
  baseURL = "http://archive.stsci.edu/cgi-bin/dss_search";

  raHours = ra / 15.;
  raH = IntegerPart[raHours];
  raM = IntegerPart[(raHours - raH) * 60];
  raS = ((raHours - raH) * 60 - raM) * 60;

  decSign = If[dec >= 0, "+", "-"];
  decAbs  = Abs[dec];
  decD = IntegerPart[decAbs];
  decM = IntegerPart[(decAbs - decD) * 60];
  decS = ((decAbs - decD) * 60 - decM) * 60;

  params = <|
    "ra"         -> StringJoin[IntegerString[raH, 10, 2], " ",
                                IntegerString[raM, 10, 2], " ",
                                StringPadLeft[ToString[NumberForm[raS, {5, 2}]], 5, "0"]],
    "dec"        -> StringJoin[decSign, IntegerString[decD, 10, 2], " ",
                                IntegerString[decM, 10, 2], " ",
                                StringPadLeft[ToString[NumberForm[decS, {5, 2}]], 5, "0"]],
    "equinox"    -> "J2000",
    "height"     -> ToString[sizeArcmin],
    "width"      -> ToString[sizeArcmin],
    "format"     -> "JPEG",
    "generation" -> "DSS2"
  |>;
  url = baseURL <> "?" <> URLQueryEncode[params];

  Print["\nQuerying STScI DSS2 for: RA=", params["ra"], ", Dec=", params["dec"]];
  Print["Field of view: ", sizeArcmin, "' x ", sizeArcmin, "'"];

  img = Quiet @ TimeConstrained[Import[url, "JPEG"], 30, $Failed];
  If[ Head[img] === Image,
    Print["\[Checkmark] Successfully retrieved DSS image"]; img,
    Print["\[Times] DSS query failed"]; $Failed
  ]
];

displayAndAnalyze[img_Image, title_String, savePath_String] := Module[
  {gray, flat, dims, row, combined, hist, histPath},
  dims = ImageDimensions[img];
  Print["\n--- Image Analysis ---"];
  Print["Image shape: ", Reverse[dims], " (plus channels ", ImageChannels[img], ")"];
  Print["Data type: ", ImageType[img]];

  gray = ColorConvert[img, "Grayscale"];
  flat = Flatten @ ImageData[gray];
  Print["Value range: [", Min[flat], ", ", Max[flat], "]"];

  combined = GraphicsGrid[{{
    Labeled[img, "Original", Top],
    Labeled[gray, "Grayscale", Top]
  }}, ImageSize -> 900];
  Export[savePath, Rasterize[combined, ImageResolution -> 150]];
  Print["\[Checkmark] Saved to: ", savePath];

  Print["\n--- Intensity Statistics (grayscale) ---"];
  Print["Mean:    ", N@Mean[flat]];
  Print["Std Dev: ", N@StandardDeviation[flat]];
  Print["Median:  ", N@Median[flat]];

  hist = Histogram[flat, 100,
    ChartStyle -> RGBColor[0.27, 0.51, 0.71],
    Frame -> True,
    FrameLabel -> {"Pixel Intensity", "Count"},
    PlotLabel -> title <> " - Intensity Distribution",
    ImageSize -> 900
  ];
  histPath = StringReplace[savePath, ".png" -> "_histogram.png"];
  Export[histPath, Rasterize[hist, ImageResolution -> 150]];
  Print["\[Checkmark] Histogram saved to: ", histPath];
];

horseheadAlternativeMain[outDir_String:"."] := Module[
  {sdss, dss},
  Print[StringJoin[ConstantArray["=", 70]]];
  Print["Astronomical Image Query: Horsehead Nebula (Barnard 33)"];
  Print[StringJoin[ConstantArray["=", 70]]];
  Print["\nCoordinates:"];
  Print["  RA:  85.2479\[Degree] = 5h 40m 59.5s"];
  Print["  Dec: -2.4583\[Degree] = -2\[Degree] 27' 30\""];
  Print["\nNote: The Horsehead Nebula is a dark nebula - an absorption feature"];
  Print["      visible against the emission nebula IC 434 in the Orion constellation."];

  Print["\n", StringJoin[ConstantArray["=", 70]]];
  Print["Method 1: SDSS DR18 (may not cover this region well)"];
  Print[StringJoin[ConstantArray["=", 70]]];
  sdss = querySDSSImage[HORSEHEAD$RA, HORSEHEAD$DEC,
            "Scale" -> 0.8, "Width" -> 800, "Height" -> 800];
  If[ Head[sdss] === Image,
    displayAndAnalyze[sdss, "Horsehead Nebula Region (SDSS)",
      FileNameJoin[{outDir, "horsehead_sdss.png"}]]
  ];

  Print["\n", StringJoin[ConstantArray["=", 70]]];
  Print["Method 2: STScI Digitized Sky Survey"];
  Print[StringJoin[ConstantArray["=", 70]]];
  dss = queryDSSViaAlternate[HORSEHEAD$RA, HORSEHEAD$DEC, 30];
  If[ Head[dss] === Image,
    displayAndAnalyze[dss, "Horsehead Nebula (DSS2)",
      FileNameJoin[{outDir, "horsehead_dss.png"}]]
  ];

  Print["\n", StringJoin[ConstantArray["=", 70]]];
  Print["Query Complete!"];
  Print[StringJoin[ConstantArray["=", 70]]];

  If[ Head[sdss] =!= Image && Head[dss] =!= Image,
    Print["\n\[WarningSign] Network restrictions prevented both queries."];
    Print["The code structure is correct - would work with proper network access."];
    Print["\nAlternative: Download FITS files manually from:"];
    Print["  - https://skyview.gsfc.nasa.gov/current/cgi/titlepage.pl"];
    Print["  - https://archive.stsci.edu/dss/"]
  ];
];

End[];
EndPackage[];

If[ !$Notebooks, HorseheadAlternative`horseheadAlternativeMain[DirectoryName[$InputFileName]] ]
