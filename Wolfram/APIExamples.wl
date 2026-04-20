(* ::Package:: *)

(* ::Title:: *)
(*APIExamples.wl*)

(* ::Text:: *)
(*Wolfram Language translation of api_examples.json.*)
(*Stores the SkyView and Las Cumbres Observatory (LCO) request templates*)
(*as a native Wolfram Association and exports them to JSON.*)


BeginPackage["APIExamples`"];

apiExamples::usage =
  "apiExamples is an Association describing example requests for the NASA SkyView \
and Las Cumbres Observatory REST APIs. Use Export[file, apiExamples, \"JSON\"] \
to serialise.";

Begin["`Private`"];

apiExamples = <|
  "SkyView" -> <|
    "url" -> "https://skyview.gsfc.nasa.gov/cgi-bin/images",
    "method" -> "GET",
    "params" -> <|
      "Position" -> "85.2479,-2.4583",
      "Survey" -> "DSS2 Red",
      "Pixels" -> "1000,1000",
      "Size" -> 0.5,
      "Return" -> "FITS"
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
      "requests" -> {
        <|
          "target" -> <|
            "name" -> "Horsehead Nebula",
            "type" -> "ICRS",
            "ra" -> 85.2479,
            "dec" -> -2.4583
          |>,
          "configurations" -> {
            <|
              "type" -> "EXPOSE",
              "instrument_type" -> "1M0-SCICAM-SINISTRO",
              "exposure_time" -> 300,
              "exposure_count" -> 3,
              "filter" -> "rp"
            |>
          },
          "windows" -> {
            <|
              "start" -> "2025-11-17T00:00:00Z",
              "end"   -> "2025-11-18T00:00:00Z"
            |>
          },
          "location" -> <|"telescope_class" -> "1m0"|>
        |>
      }
    |>,
    "wolfram" -> "
request = HTTPRequest[
  \"https://observe.lco.global/api/requestgroups/\",
  <|
    \"Method\"  -> \"POST\",
    \"Headers\" -> {\"Authorization\" -> \"Token YOUR_API_TOKEN\"},
    \"ContentType\" -> \"application/json\",
    \"Body\"    -> ExportString[observationRequest, \"JSON\"]
  |>];
response = URLRead[request];
"
  |>
|>;

exportAPIExamples[file_String] := Export[file, apiExamples, "JSON"];

End[];
EndPackage[];

(* ::Section:: *)
(*Script mode: write api_examples.json next to this file*)

If[ !$Notebooks,
  Module[{outFile},
    outFile = FileNameJoin[{DirectoryName[$InputFileName], "api_examples.json"}];
    Export[outFile, APIExamples`apiExamples, "JSON"];
    Print["Wrote ", outFile];
  ]
]
