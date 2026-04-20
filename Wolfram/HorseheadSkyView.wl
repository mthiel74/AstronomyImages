(* ::Package:: *)

(* ::Title:: *)
(*HorseheadSkyView.wl*)

(* ::Text:: *)
(*Wolfram Language translation of horsehead_skyview.py.*)
(*Queries the NASA SkyView Virtual Observatory for the Horsehead Nebula*)
(*across several surveys, saves the FITS and a PNG rendering of each.*)


BeginPackage["HorseheadSkyView`"];

queryskyView::usage =
  "querySkyView[ra, dec, opts] queries NASA SkyView for a FITS cutout at the given \
equatorial coordinates. Options: \"Survey\" (\"DSS2 Red\"), \"Pixels\" (800), \"FOV\" (0.5).";

displayFITSImage::usage =
  "displayFITSImage[fitsFile, title] returns a rescaled Image constructed from the \
FITS primary HDU with percentile clipping. Prints basic header statistics.";

saveFITS::usage =
  "saveFITS[bytes, file] writes the raw FITS byte string to disk.";

horseheadSkyViewMain::usage =
  "horseheadSkyViewMain[] replicates the Python main(): query DSS2 Red and 2MASS-J, \
save FITS and PNG images for each.";

Begin["`Private`"];

(* Horsehead Nebula coordinates (Barnard 33 in Orion) *)
(* RA 05h 40m 59s = 85.2479 deg, Dec -02\[Degree] 27' 30" = -2.4583 deg *)
HORSEHEAD$RA  = 85.2479;
HORSEHEAD$DEC = -2.4583;

Options[querySkyView] = {
  "Survey"  -> "DSS2 Red",
  "Pixels"  -> 800,
  "FOV"     -> 0.5,
  "Timeout" -> 60
};

querySkyView[ra_?NumericQ, dec_?NumericQ, OptionsPattern[]] := Module[
  {survey, pixels, fov, timeout, baseURL, params, url, bytes},
  survey  = OptionValue["Survey"];
  pixels  = OptionValue["Pixels"];
  fov     = OptionValue["FOV"];
  timeout = OptionValue["Timeout"];

  baseURL = "https://skyview.gsfc.nasa.gov/cgi-bin/images";
  params = <|
    "Position" -> ToString[ra] <> "," <> ToString[dec],
    "Survey"   -> survey,
    "Pixels"   -> ToString[pixels] <> "," <> ToString[pixels],
    "Size"     -> ToString[fov],
    "Return"   -> "FITS",
    "Scaling"  -> "Linear"
  |>;
  url = baseURL <> "?" <> URLQueryEncode[params];

  Print["Querying SkyView for coordinates: RA=", ra, "\[Degree], Dec=", dec, "\[Degree]"];
  Print["Survey: ", survey, ", FOV: ", fov, "\[Degree], Resolution: ",
        pixels, "x", pixels, " pixels"];

  bytes = TimeConstrained[
    URLExecute[url, "ByteArray"],
    timeout,
    $Failed
  ];
  If[ bytes === $Failed || bytes === {} || ByteArray[] === bytes,
    Throw["SkyView query failed", querySkyView]
  ];
  Print["\[Checkmark] Successfully retrieved data (", ByteCount[bytes], " bytes)"];
  bytes
];

saveFITS[bytes_ByteArray, filename_String] := (
  BinaryWrite[filename, bytes];
  Close[filename];
  Print["\[Checkmark] FITS file saved to: ", filename];
  filename
);

(* Display a FITS file (path) with percentile clipping. Returns an Image. *)
displayFITSImage[fitsFile_String, title_String:"Astronomical Image"] := Module[
  {data, meta, dims, vmin, vmax, scaled, img},
  data = Quiet @ Import[fitsFile, {"FITS", "Data"}];
  If[ Head[data] === List && Depth[data] > 2, data = First[data] ];
  meta = Quiet @ Import[fitsFile, "Metadata"];

  If[ !MatrixQ[data, NumericQ],
    Print["Could not decode FITS numeric data from ", fitsFile];
    Return[$Failed]
  ];

  dims = Dimensions[data];
  Print["\n--- FITS Header Information ---"];
  Print["Image shape: ", dims];
  Print["Min value:  ", ScientificForm[N@Min[DeleteCases[Flatten[data], _?(Not@*NumericQ)]]]];
  Print["Max value:  ", ScientificForm[N@Max[DeleteCases[Flatten[data], _?(Not@*NumericQ)]]]];
  Print["Mean value: ", ScientificForm[N@Mean[DeleteCases[Flatten[data], _?(Not@*NumericQ)]]]];

  {vmin, vmax} = Quantile[DeleteCases[Flatten[data], _?(Not@*NumericQ)], {0.01, 0.995}];
  scaled = Clip[(data - vmin) / (vmax - vmin + 10.^-12), {0, 1}];
  img = Image[Reverse[scaled], "Real", ColorFunction -> GrayLevel];

  Labeled[img, Style[title, Bold, 14]]
];

downloadAndPlot[ra_, dec_, survey_, description_, fov_, pixels_, outDir_] := Module[
  {bytes, safe, fitsPath, pngPath, img},
  Print[StringJoin[ConstantArray["=", 60]]];
  Print["Querying: ", description];
  Print[StringJoin[ConstantArray["=", 60]]];

  bytes = Catch[querySkyView[ra, dec,
    "Survey" -> survey, "Pixels" -> pixels, "FOV" -> fov], querySkyView];
  If[ !ByteArrayQ[bytes],
    Print["\[Times] Error querying ", survey, ": ", bytes];
    Return[$Failed]
  ];

  safe = StringReplace[survey, {" " -> "_", "-" -> "_"}];
  fitsPath = FileNameJoin[{outDir, "horsehead_" <> safe <> ".fits"}];
  pngPath  = FileNameJoin[{outDir, "horsehead_" <> safe <> ".png"}];

  saveFITS[bytes, fitsPath];
  img = displayFITSImage[fitsPath, "Horsehead Nebula - " <> description];
  If[ Head[img] =!= $Failed,
    Export[pngPath, img];
    Print["\[Checkmark] Plot saved to: ", pngPath];
  ];
  img
];

horseheadSkyViewMain[outDir_String:"."] := Module[
  {surveys},
  Print[StringJoin[ConstantArray["=", 60]]];
  Print["NASA SkyView Query: Horsehead Nebula (Barnard 33)"];
  Print[StringJoin[ConstantArray["=", 60]]];

  surveys = {
    {"DSS2 Red", "Digitized Sky Survey (Red)", 0.5},
    {"2MASS-J",  "2MASS Near-Infrared (J-band)", 0.5}
  };

  Scan[
    Function[entry,
      downloadAndPlot[HORSEHEAD$RA, HORSEHEAD$DEC,
        entry[[1]], entry[[2]], entry[[3]], 1000, outDir]
    ],
    surveys
  ];

  Print["\n", StringJoin[ConstantArray["=", 60]]];
  Print["Analysis complete!"];
  Print[StringJoin[ConstantArray["=", 60]]];
  Print["\n--- Statistical Analysis ---"];
  Print["The Horsehead Nebula is a dark nebula (absorption feature)"];
  Print["Expected: Lower intensity in nebula region vs. surrounding emission"];
];

End[];
EndPackage[];

(* Script mode entry point *)
If[ !$Notebooks, HorseheadSkyView`horseheadSkyViewMain[DirectoryName[$InputFileName]] ]
