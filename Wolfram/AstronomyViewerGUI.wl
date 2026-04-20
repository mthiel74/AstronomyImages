(* ::Package:: *)

(* ::Title:: *)
(*AstronomyViewerGUI.wl*)

(* ::Text:: *)
(*Wolfram Language translation of astronomy_viewer_gui.py.*)
(*Replaces the Tk-based Python UI with a DynamicModule so you can pick a*)
(*catalogue object or enter coordinates, download a FITS cutout from NASA*)
(*SkyView, and inspect/save the image directly in a Wolfram notebook.*)


BeginPackage["AstronomyViewerGUI`"];

catalog::usage =
  "catalog is an Association of named astronomical objects with RA, Dec, and angular size.";

launchAstronomyViewer::usage =
  "launchAstronomyViewer[] opens an interactive DynamicModule in the front end. \
In script mode it returns the expression so it can be evaluated manually.";

Begin["`Private`"];

catalog = <|
  "Horsehead Nebula"       -> <|"ra" -> 85.2479,  "dec" -> -2.4583,  "size" -> 0.3|>,
  "Orion Nebula (M42)"     -> <|"ra" -> 83.8221,  "dec" -> -5.3911,  "size" -> 1.0|>,
  "Andromeda Galaxy (M31)" -> <|"ra" -> 10.6847,  "dec" -> 41.2687,  "size" -> 3.0|>,
  "Crab Nebula (M1)"       -> <|"ra" -> 83.6333,  "dec" -> 22.0145,  "size" -> 0.3|>,
  "Ring Nebula (M57)"      -> <|"ra" -> 283.3963, "dec" -> 33.0295,  "size" -> 0.1|>,
  "Whirlpool Galaxy (M51)" -> <|"ra" -> 202.4696, "dec" -> 47.1952,  "size" -> 0.3|>,
  "Eagle Nebula (M16)"     -> <|"ra" -> 274.7000, "dec" -> -13.8167, "size" -> 0.5|>,
  "Lagoon Nebula (M8)"     -> <|"ra" -> 270.9167, "dec" -> -24.3833, "size" -> 1.0|>,
  "Pleiades (M45)"         -> <|"ra" -> 56.75,    "dec" -> 24.1167,  "size" -> 2.0|>,
  "Dumbbell Nebula (M27)"  -> <|"ra" -> 299.9013, "dec" -> 22.7211,  "size" -> 0.3|>
|>;

(* Pull a FITS byte array from SkyView *)
skyViewQuery[ra_, dec_, survey_, pixels_, fov_] := Module[
  {url, params},
  url = "https://skyview.gsfc.nasa.gov/cgi-bin/images";
  params = <|
    "Position" -> ToString[ra] <> "," <> ToString[dec],
    "Survey"   -> survey,
    "Pixels"   -> ToString[pixels] <> "," <> ToString[pixels],
    "Size"     -> ToString[fov],
    "Return"   -> "FITS",
    "Scaling"  -> "Linear"
  |>;
  TimeConstrained[URLExecute[url <> "?" <> URLQueryEncode[params], "ByteArray"], 60, $Failed]
];

(* Render the numeric image data with a given ColorFunction name. *)
renderImage[data_?MatrixQ, cmapName_] := Module[
  {flat, vmin, vmax, scaled, cf},
  flat = DeleteCases[Flatten[data], _?(Not@*NumericQ)];
  {vmin, vmax} = Quantile[flat, {0.01, 0.99}];
  scaled = Clip[(data - vmin) / (vmax - vmin + 10.^-12), {0, 1}];
  cf = ColorData[cmapName];
  Image[Reverse[scaled], "Real", ColorFunction -> cf]
];

(* Plot the image with a colorbar via ArrayPlot. *)
arrayPlotWithBar[data_, cmapName_, title_] := Module[
  {flat, vmin, vmax},
  flat = DeleteCases[Flatten[data], _?(Not@*NumericQ)];
  {vmin, vmax} = Quantile[flat, {0.01, 0.99}];
  ArrayPlot[Reverse[data],
    ColorFunction -> ColorData[cmapName],
    ColorFunctionScaling -> True,
    PlotRange -> {vmin, vmax},
    Frame -> True,
    FrameLabel -> {"X (pixels)", "Y (pixels)"},
    PlotLabel -> Style[title, Bold, 12],
    PlotLegends -> Automatic,
    ImageSize -> 620
  ]
];

launchAstronomyViewer[] := DynamicModule[
  {objectName = "Horsehead Nebula",
   raStr = "85.2479", decStr = "-2.4583", fovStr = "0.5",
   survey = "DSS2 Red", resolution = "800", cmap = "GrayTones",
   status = "Ready",
   currentBytes = None, currentData = None, info = "No image loaded.",
   objectNames, surveys, resolutions, cmaps, display},

  objectNames = Keys[catalog];
  surveys     = {"DSS2 Red", "DSS2 Blue", "2MASS-J", "2MASS-H", "WISE 3.4"};
  resolutions = {"500", "800", "1000", "1500", "2000"};
  cmaps       = {"GrayTones", "Rainbow", "SunsetColors", "TemperatureMap",
                 "DeepSeaColors", "AuroraColors"};

  display = Function[{},
    If[ currentData =!= None,
      arrayPlotWithBar[currentData, cmap,
        survey <> " - RA: " <> raStr <> "\[Degree], Dec: " <> decStr <> "\[Degree]"],
      Graphics[Text[Style["Select a catalog object or enter\ncoordinates, then Download Image",
        14], {0, 0}], ImageSize -> 620]
    ]
  ];

  Grid[{
    { (* Controls *)
      Panel[Column[{
        Style["Controls", Bold, 14],
        Row[{"Catalog:",
          PopupMenu[Dynamic[objectName], objectNames, FieldSize -> 25]
        }, Spacer[4]],
        Button["Load Object",
          Module[{obj},
            obj = catalog[objectName];
            raStr  = ToString[obj["ra"]];
            decStr = ToString[obj["dec"]];
            fovStr = ToString[obj["size"]];
            info   = "Loaded: " <> objectName;
          ], ImageSize -> 200],
        Delimiter,
        Style["Manual Coordinates", Bold],
        Row[{"RA (\[Degree]):", InputField[Dynamic[raStr], String, FieldSize -> 12]}, Spacer[4]],
        Row[{"Dec (\[Degree]):", InputField[Dynamic[decStr], String, FieldSize -> 12]}, Spacer[4]],
        Row[{"FOV (\[Degree]):", InputField[Dynamic[fovStr], String, FieldSize -> 12]}, Spacer[4]],
        Row[{"Survey:", PopupMenu[Dynamic[survey], surveys]}, Spacer[4]],
        Row[{"Resolution:", PopupMenu[Dynamic[resolution], resolutions]}, Spacer[4]],
        Button[Style["Download Image", Bold, Blue],
          Module[{bytes, tmp, data},
            status = "Downloading...";
            bytes = skyViewQuery[ToExpression[raStr], ToExpression[decStr],
                                 survey, ToExpression[resolution],
                                 ToExpression[fovStr]];
            If[ ByteArrayQ[bytes] && ByteCount[bytes] > 1000,
              currentBytes = bytes;
              tmp = CreateFile[];
              BinaryWrite[tmp, bytes]; Close[tmp];
              data = Quiet @ Import[tmp, {"FITS", "Data"}];
              If[ Head[data] === List && Depth[data] > 2, data = First[data]];
              If[ MatrixQ[data, NumericQ],
                currentData = data;
                info = "Survey: " <> survey <>
                       "\nCoordinates: RA=" <> raStr <> "\[Degree], Dec=" <> decStr <> "\[Degree]" <>
                       "\nImage size: " <> ToString[Dimensions[data]] <>
                       "\nIntensity range: [" <>
                         ToString[NumberForm[Min[Flatten[data]], {5, 2}]] <> ", " <>
                         ToString[NumberForm[Max[Flatten[data]], {5, 2}]] <> "]";
                status = "\[Checkmark] Download complete",
                status = "\[Times] Could not decode FITS"
              ],
              status = "\[Times] Download failed"
            ];
          ], ImageSize -> 200],
        Delimiter,
        Style["Display", Bold],
        Row[{"Color Map:", PopupMenu[Dynamic[cmap], cmaps]}, Spacer[4]],
        Button["Save FITS",
          If[ currentBytes =!= None,
            Module[{path}, path = SystemDialogInput["FileSave", {Directory[], "*.fits"}];
              If[ StringQ[path], BinaryWrite[path, currentBytes]; Close[path];
                info = "Saved FITS to " <> path]
            ],
            info = "No image downloaded yet."
          ], ImageSize -> 200],
        Button["Save PNG",
          If[ currentData =!= None,
            Module[{path}, path = SystemDialogInput["FileSave", {Directory[], "*.png"}];
              If[ StringQ[path],
                Export[path, Rasterize[display[], ImageResolution -> 150]];
                info = "Saved PNG to " <> path]
            ],
            info = "No image downloaded yet."
          ], ImageSize -> 200],
        Dynamic[Style["Status: " <> status, Italic]]
      }, Spacings -> 1], ImageSize -> 260],

      (* Display *)
      Panel[Column[{
        Style["Image Information", Bold, 14],
        Dynamic[info],
        Dynamic[display[]]
      }, Spacings -> 1], ImageSize -> {680, 720}]
    }
  }, Alignment -> Top, Dividers -> Center]
];

End[];
EndPackage[];

If[ !$Notebooks,
  Print["AstronomyViewerGUI requires a Wolfram front end (notebook). "];
  Print["In a notebook, evaluate: Get[\"",
        FileNameJoin[{DirectoryName[$InputFileName], "AstronomyViewerGUI.wl"}], "\"]"];
  Print["then AstronomyViewerGUI`launchAstronomyViewer[]"]
]
