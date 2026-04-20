(* ::Package:: *)

(* ::Title:: *)
(*HorseheadSkyViewFixed.wl*)

(* ::Text:: *)
(*Wolfram Language translation of horsehead_skyview_fixed.py.*)
(*Same as HorseheadSkyView.wl but writes FITS/PNG into the current working*)
(*directory (no /mnt/user-data/outputs prefix).*)


Needs["HorseheadSkyView`", FileNameJoin[{DirectoryName[$InputFileName], "HorseheadSkyView.wl"}]];

BeginPackage["HorseheadSkyViewFixed`", {"HorseheadSkyView`"}];

horseheadSkyViewFixedMain::usage =
  "horseheadSkyViewFixedMain[] queries SkyView for DSS2 Red and 2MASS-J cutouts of \
the Horsehead Nebula and saves FITS + PNG alongside this script.";

Begin["`Private`"];

HORSEHEAD$RA  = 85.2479;
HORSEHEAD$DEC = -2.4583;

horseheadSkyViewFixedMain[] := Module[
  {surveys, dir},
  dir = If[ StringQ[$InputFileName] && $InputFileName =!= "",
           DirectoryName[$InputFileName],
           Directory[] ];

  Print[StringJoin[ConstantArray["=", 60]]];
  Print["NASA SkyView Query: Horsehead Nebula (Barnard 33)"];
  Print[StringJoin[ConstantArray["=", 60]]];

  surveys = {
    {"DSS2 Red", "Digitized Sky Survey (Red)", 0.5},
    {"2MASS-J",  "2MASS Near-Infrared (J-band)", 0.5}
  };

  Scan[
    Function[entry, Module[{survey, description, fov, safe, bytes, fitsPath, pngPath, img},
      survey = entry[[1]]; description = entry[[2]]; fov = entry[[3]];
      Print["\n", StringJoin[ConstantArray["=", 60]]];
      Print["Querying: ", description];
      Print[StringJoin[ConstantArray["=", 60]]];

      bytes = Catch[
        HorseheadSkyView`querySkyView[HORSEHEAD$RA, HORSEHEAD$DEC,
          "Survey" -> survey, "Pixels" -> 1000, "FOV" -> fov],
        HorseheadSkyView`querySkyView];

      If[ !ByteArrayQ[bytes],
        Print["\[Times] Error querying ", survey, ": ", bytes];
        Return[Null]
      ];

      safe     = StringReplace[survey, {" " -> "_", "-" -> "_"}];
      fitsPath = FileNameJoin[{dir, "horsehead_" <> safe <> ".fits"}];
      pngPath  = FileNameJoin[{dir, "horsehead_" <> safe <> ".png"}];

      HorseheadSkyView`saveFITS[bytes, fitsPath];
      img = HorseheadSkyView`displayFITSImage[fitsPath,
              "Horsehead Nebula - " <> description];
      If[ Head[img] =!= $Failed,
        Export[pngPath, img];
        Print["\[Checkmark] Plot saved to: ", pngPath]
      ];
    ]],
    surveys
  ];

  Print["\n", StringJoin[ConstantArray["=", 60]]];
  Print["Analysis complete!"];
  Print[StringJoin[ConstantArray["=", 60]]];
];

End[];
EndPackage[];

If[ !$Notebooks, HorseheadSkyViewFixed`horseheadSkyViewFixedMain[] ]
