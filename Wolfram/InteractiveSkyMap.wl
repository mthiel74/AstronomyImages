(* ::Package:: *)

(* ::Title:: *)
(*InteractiveSkyMap.wl*)

(* ::Text:: *)
(*Wolfram Language translation of interactive_sky_map.py.*)
(*A DynamicModule with a clickable all-sky map, a controls panel and a*)
(*telescope view panel. Click anywhere on the sky map, adjust the FOV /*)
(*survey / resolution, and download a FITS cutout from NASA SkyView.*)


BeginPackage["InteractiveSkyMap`"];

constellations::usage =
  "constellations is an Association of (approximate) constellation centers.";

famousObjects::usage =
  "famousObjects is an Association of selected deep-sky targets (RA, Dec, label).";

launchInteractiveSkyMap::usage =
  "launchInteractiveSkyMap[] returns a DynamicModule expression providing a clickable \
all-sky map with a SkyView FITS downloader and viewer.";

Begin["`Private`"];

constellations = <|
  "Orion"       -> <|"ra" -> 83.8,  "dec" -> -5.4,  "size" -> 30|>,
  "Ursa Major"  -> <|"ra" -> 165.0, "dec" -> 55.0,  "size" -> 25|>,
  "Cassiopeia"  -> <|"ra" -> 15.0,  "dec" -> 60.0,  "size" -> 20|>,
  "Andromeda"   -> <|"ra" -> 10.7,  "dec" -> 41.3,  "size" -> 20|>,
  "Cygnus"      -> <|"ra" -> 305.0, "dec" -> 40.0,  "size" -> 25|>,
  "Sagittarius" -> <|"ra" -> 283.0, "dec" -> -25.0, "size" -> 30|>,
  "Taurus"      -> <|"ra" -> 68.0,  "dec" -> 15.0,  "size" -> 25|>,
  "Leo"         -> <|"ra" -> 152.0, "dec" -> 15.0,  "size" -> 25|>,
  "Aquila"      -> <|"ra" -> 297.0, "dec" -> 10.0,  "size" -> 20|>,
  "Centaurus"   -> <|"ra" -> 192.0, "dec" -> -43.0, "size" -> 30|>
|>;

famousObjects = <|
  "Horsehead" -> <|"ra" -> 85.2479,  "dec" -> -2.4583, "label" -> "Horsehead"|>,
  "M31"       -> <|"ra" -> 10.6847,  "dec" -> 41.2687, "label" -> "M31 (Andromeda)"|>,
  "M42"       -> <|"ra" -> 83.8221,  "dec" -> -5.3911, "label" -> "M42 (Orion)"|>,
  "M51"       -> <|"ra" -> 202.4696, "dec" -> 47.1952, "label" -> "M51 (Whirlpool)"|>,
  "M57"       -> <|"ra" -> 283.3963, "dec" -> 33.0295, "label" -> "M57 (Ring)"|>,
  "Crab"      -> <|"ra" -> 83.6333,  "dec" -> 22.0145, "label" -> "M1 (Crab)"|>
|>;

(* Convert RA (deg) and Dec (deg) to sexagesimal strings. *)
raToHMS[ra_?NumericQ] := Module[{raH, h, m, s},
  raH = Mod[ra, 360] / 15.;
  h = IntegerPart[raH];
  m = IntegerPart[(raH - h) * 60];
  s = ((raH - h) * 60 - m) * 60;
  IntegerString[h, 10, 2] <> ":" <> IntegerString[m, 10, 2] <> ":" <>
    ToString[NumberForm[s, {4, 1}]]
];

decToDMS[dec_?NumericQ] := Module[{sign, a, d, m, s},
  sign = If[dec >= 0, "+", "-"];
  a = Abs[dec];
  d = IntegerPart[a];
  m = IntegerPart[(a - d) * 60];
  s = ((a - d) * 60 - m) * 60;
  sign <> IntegerString[d, 10, 2] <> ":" <> IntegerString[m, 10, 2] <> ":" <>
    ToString[NumberForm[s, {4, 1}]]
];

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
  TimeConstrained[
    Module[{resp},
      resp = URLRead[HTTPRequest[url <> "?" <> URLQueryEncode[params]]];
      If[resp === $Failed || resp["StatusCode"] =!= 200, $Failed,
        resp["BodyByteArray"]]
    ],
    60, $Failed]
];

(* Build the static all-sky map graphics, with an optional "Selected" mark. *)
skyMap[selRA_, selDec_] := Module[
  {raLine, ecliptic, mwRA, mwDec, constPatches, objPoints, selection},
  raLine   = Subdivide[0., 360., 99];
  ecliptic = {#, 23.4 Sin[# Degree]} & /@ raLine;
  mwRA     = Subdivide[0., 360., 99];
  mwDec    = -29 + 60 Sin[(# - 266) Degree] & /@ mwRA;

  constPatches = KeyValueMap[
    Function[{name, data},
      {Opacity[0.1], Cyan, Disk[{data["ra"], data["dec"]}, data["size"]/2],
       Opacity[0.7], Darker[Cyan],
       Text[Style[name, Bold, 8], {data["ra"], data["dec"]}]}
    ],
    constellations];

  objPoints = KeyValueMap[
    Function[{name, data},
      {Red, PointSize[Large], Point[{data["ra"], data["dec"]}],
       Text[Style[data["label"], Bold, 8, Red],
            {data["ra"], data["dec"] + 5}, {0, -1}]}
    ],
    famousObjects];

  selection = If[ NumericQ[selRA] && NumericQ[selDec],
    {Green, PointSize[Large], Point[{selRA, selDec}],
     Yellow, Text[Style["SELECTED", Bold, 9], {selRA, selDec - 8}, {0, 1}]},
    {}];

  Graphics[{
    (* Celestial equator *)
    {Blue, Opacity[0.3], Thickness[0.004], Line[{#, 0} & /@ raLine]},
    (* Ecliptic *)
    {Orange, Opacity[0.3], Dashed, Line[ecliptic]},
    (* Milky Way *)
    {Purple, Opacity[0.2], Thickness[0.01], Line[MapThread[List, {mwRA, mwDec}]]},
    Sequence @@ constPatches,
    Sequence @@ objPoints,
    Sequence @@ selection
  },
    Frame -> True,
    FrameLabel -> {"Right Ascension (\[Degree])", "Declination (\[Degree])"},
    PlotLabel -> Style["All-Sky Map (Click to Select)", Bold, 12],
    PlotRange -> {{360, 0}, {-90, 90}},  (* RA increases to the left *)
    AspectRatio -> 1/2,
    GridLines -> Automatic,
    GridLinesStyle -> Directive[Gray, Opacity[0.3], Dashed],
    ImageSize -> 640
  ]
];

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
    PlotLabel -> Style[title, Bold, 11],
    PlotLegends -> Automatic,
    ImageSize -> 640
  ]
];

launchInteractiveSkyMap[] := DynamicModule[
  {selRA = None, selDec = None,
   fov = "0.5", survey = "DSS2 Red", resolution = "1000", cmap = "GrayTones",
   status = "Ready - click sky map",
   info = "Click anywhere on the sky map to select coordinates,\nthen Download & View to observe the region.",
   currentBytes = None, currentData = None,
   surveys, resolutions, cmaps},

  surveys = {"DSS2 Red", "DSS2 Blue", "DSS2 IR", "2MASS-J", "2MASS-H",
             "2MASS-K", "WISE 3.4"};
  resolutions = {"500", "800", "1000", "1500", "2000"};
  cmaps = {"GrayTones", "Rainbow", "SunsetColors", "TemperatureMap",
           "DeepSeaColors", "AuroraColors"};

  Grid[{{
    (* LEFT: clickable sky map *)
    Panel[Column[{
      Style["\[CircleTimes] Interactive Sky Map - Click Anywhere!", Bold, 12],
      EventHandler[
        Dynamic[skyMap[selRA, selDec]],
        {"MouseClicked" :> (
          {selRA, selDec} = MousePosition["Graphics"];
          (* MousePosition returns $Failed if outside graphics; guard. *)
          If[ Head[selRA] === Real || IntegerQ[selRA],
            status = "Coordinates selected - ready to download";
            info = "Selected: RA=" <> ToString[NumberForm[selRA, {6, 4}]] <>
                   "\[Degree] (" <> raToHMS[selRA] <> "), Dec=" <>
                   ToString[NumberForm[selDec, {6, 4}]] <>
                   "\[Degree] (" <> decToDMS[selDec] <> ")\n" <>
                   "Click 'Download & View' to observe this region.",
            selRA = None; selDec = None
          ]
        )}
      ],
      Dynamic @ Style[
        "Selected: " <>
          If[ NumericQ[selRA],
            "RA=" <> ToString[NumberForm[selRA, {6, 4}]] <> "\[Degree], " <>
            "Dec=" <> ToString[NumberForm[selDec, {6, 4}]] <> "\[Degree]",
            "(none)"], Italic]
    }, Spacings -> 1], ImageSize -> {680, 720}],

    (* MIDDLE: controls *)
    Panel[Column[{
      Style["Controls", Bold, 12],
      Row[{"FOV (\[Degree]):", InputField[Dynamic[fov], String, FieldSize -> 8]}, Spacer[4]],
      Row[{
        Button["0.1\[Degree]", fov = "0.1", ImageSize -> 50],
        Button["0.5\[Degree]", fov = "0.5", ImageSize -> 50],
        Button["1.0\[Degree]", fov = "1.0", ImageSize -> 50],
        Button["2.0\[Degree]", fov = "2.0", ImageSize -> 50]
      }],
      Row[{"Survey:",     PopupMenu[Dynamic[survey], surveys]}, Spacer[4]],
      Row[{"Resolution:", PopupMenu[Dynamic[resolution], resolutions]}, Spacer[4]],
      Button[Style["\[DoubleVerticalBar] Download & View", Bold, Blue],
        If[ NumericQ[selRA] && NumericQ[selDec],
          Module[{bytes, tmp, data},
            status = "Downloading from SkyView...";
            bytes = skyViewQuery[selRA, selDec, survey,
              ToExpression[resolution], ToExpression[fov]];
            If[ ByteArrayQ[bytes] && ByteCount[bytes] > 1000,
              currentBytes = bytes;
              tmp = CreateFile[];
              BinaryWrite[tmp, bytes]; Close[tmp];
              data = Quiet @ Import[tmp, {"FITS", "Data"}];
              data = Which[
                AssociationQ[data],                   First[Values[data]],
                ListQ[data] && Depth[data] >= 4,      First[data],
                True,                                 data
              ];
              If[ MatrixQ[data, NumericQ],
                currentData = data;
                info = "Survey: " <> survey <> " | FOV: " <> fov <>
                       "\[Degree] | Resolution: " <> ToString[Dimensions[data]] <>
                       "\nPosition: RA=" <> ToString[NumberForm[selRA, {6, 4}]] <>
                       "\[Degree] (" <> raToHMS[selRA] <> "), Dec=" <>
                       ToString[NumberForm[selDec, {6, 4}]] <> "\[Degree] (" <>
                       decToDMS[selDec] <> ")\nIntensity: [" <>
                       ToString[NumberForm[Min[Flatten[data]], {5, 2}]] <> ", " <>
                       ToString[NumberForm[Max[Flatten[data]], {5, 2}]] <> "]";
                status = "\[Checkmark] Download complete!",
                status = "\[Times] Could not decode FITS"
              ],
              status = "\[Times] Download failed"
            ]
          ],
          status = "Select coordinates on the sky map first."
        ], ImageSize -> 220],
      Delimiter,
      Style["Display", Bold],
      Row[{"Color map:", PopupMenu[Dynamic[cmap], cmaps]}, Spacer[4]],
      Button["Save FITS",
        If[ currentBytes =!= None,
          Module[{path}, path = SystemDialogInput["FileSave",
              {Directory[],
               "sky_RA" <> ToString[NumberForm[selRA, {5, 2}]] <> "_Dec" <>
                 ToString[NumberForm[selDec, {5, 2}]] <> ".fits"}];
            If[ StringQ[path], BinaryWrite[path, currentBytes]; Close[path];
              status = "Saved FITS"]
          ],
          status = "No image downloaded yet."
        ], ImageSize -> 220],
      Button["Save PNG",
        If[ currentData =!= None,
          Module[{path, fig},
            path = SystemDialogInput["FileSave",
              {Directory[],
               "sky_RA" <> ToString[NumberForm[selRA, {5, 2}]] <> "_Dec" <>
                 ToString[NumberForm[selDec, {5, 2}]] <> ".png"}];
            If[ StringQ[path],
              fig = arrayPlotWithBar[currentData, cmap,
                survey <> " - FOV: " <> fov <> "\[Degree]"];
              Export[path, Rasterize[fig, ImageResolution -> 150]];
              status = "Saved PNG"]
          ],
          status = "No image downloaded yet."
        ], ImageSize -> 220],
      Delimiter,
      Dynamic[Style["Status: " <> status, Italic]]
    }, Spacings -> 1], ImageSize -> 280],

    (* RIGHT: telescope view *)
    Panel[Column[{
      Style["\[Telephone] Telescope View", Bold, 12],
      Dynamic[info],
      Dynamic[
        If[ currentData =!= None,
          arrayPlotWithBar[currentData, cmap,
            survey <> " - FOV: " <> fov <> "\[Degree]\n" <>
            "RA: " <> raToHMS[selRA] <> " (" <> ToString[NumberForm[selRA, {6, 4}]] <>
            "\[Degree]), Dec: " <> decToDMS[selDec] <> " (" <>
            ToString[NumberForm[selDec, {6, 4}]] <> "\[Degree])"],
          Graphics[Text[Style[
            "Telescope View\n\nClick on the sky map to select\na region, then Download & View.",
            14], {0, 0}], ImageSize -> 640]
        ]
      ]
    }, Spacings -> 1], ImageSize -> {680, 720}]
  }}, Alignment -> Top, Dividers -> Center]
];

End[];
EndPackage[];

If[ !$Notebooks,
  Print["InteractiveSkyMap requires a Wolfram front end (notebook)."];
  Print["Evaluate: Get[\"",
        FileNameJoin[{DirectoryName[$InputFileName], "InteractiveSkyMap.wl"}], "\"]"];
  Print["then InteractiveSkyMap`launchInteractiveSkyMap[]"]
]
