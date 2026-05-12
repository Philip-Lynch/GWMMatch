(* =========================================================== *)
(*  02-PSDs.wls                                                *)
(*  GWMMatch — exploring built-in detector PSD models         *)
(*                                                             *)
(*  Run:  wolframscript -file 02-PSDs.wls                     *)
(*  Or:   Get["path/to/notebooks/02-PSDs.wls"]               *)
(*                                                             *)
(*  Plots are saved to 02-PSDs-output/ as PDF files           *)
(*  when run from the command line, and displayed inline       *)
(*  when evaluated in a Mathematica notebook.                  *)
(* =========================================================== *)

<< GWMMatch`

outputDir = FileNameJoin[{DirectoryName[$InputFileName], "02-PSDs-output"}];
If[!DirectoryQ[outputDir], CreateDirectory[outputDir]];

(* ----------------------------------------------------------- *)
(* Grid parameters shared by all plots                         *)
(* ----------------------------------------------------------- *)

deltaF = 0.5;   (* Hz bin spacing *)
fLow   = 1.0;   (* Hz — low enough to cover LISA and ET *)
nFreq  = 8193;  (* covers up to 4096 Hz at 0.5 Hz spacing *)

freqs = Table[(k - 1) * deltaF, {k, 2, nFreq}];  (* skip DC bin *)

(* Helper: extract in-band ASD values *)
psdToASD[psd_] := Sqrt[Select[
  Transpose[{freqs, Rest[psd]}],
  (#[[2]] =!= Infinity && #[[2]] > 0) &
]];

(* ----------------------------------------------------------- *)
(* 1. All available models                                     *)
(* ----------------------------------------------------------- *)

allModels = ListDetectorPSDs[];
Print["Available PSD models (", Length[allModels], "):"];
Print[Column[allModels]];

(* ----------------------------------------------------------- *)
(* 2. aLIGO: analytical vs tabulated comparison               *)
(* ----------------------------------------------------------- *)

Print["\n-- aLIGO: analytical vs GWINC tabulated --"];

psdAnalytical  = DetectorPSD["aLIGOZeroDetHighPower",        nFreq, deltaF, fLow];
psdGWINC       = DetectorPSD["aLIGOZeroDetHighPowerGWINC",   nFreq, deltaF, fLow];

asdAnalytical  = psdToASD[psdAnalytical];
asdGWINC       = psdToASD[psdGWINC];

plotALIGO = LogLogPlot[
  {Sqrt[PSDaLIGOZeroDetHighPower[f]],
   Sqrt[GWMMatch`Private`loadTabulatedASD["LIGO-T0900288-v3-ZERO_DET_high_P.txt"][f]]},
  {f, 10, 2000},
  PlotStyle -> {{Blue, Thick}, {Red, Dashed, Thick}},
  PlotLegends -> Placed[
    LineLegend[{"Analytical (LALSim)", "GWINC tabulated"},
               LegendMarkerSize -> 30],
    {0.7, 0.25}],
  PlotLabel -> "aLIGO Zero-Detuned High Power: Analytical vs Tabulated",
  AxesLabel -> {"f [Hz]", "ASD [1/\[Sqrt]Hz]"},
  PlotRange -> {{10, 2000}, {1*^-24, 1*^-21}},
  GridLines -> Automatic, GridLinesStyle -> LightGray,
  ImageSize -> 600
];

Print["aLIGO comparison plot generated."];
Export[FileNameJoin[{outputDir, "aLIGO_comparison.pdf"}], plotALIGO];

(* ----------------------------------------------------------- *)
(* 3. Next-generation ground-based: ET and CE variants        *)
(* ----------------------------------------------------------- *)

Print["-- Next-generation ground-based PSDs --"];

ngModels = {"EinsteinTelescopeP1600143",
            "CosmicExplorerP1600143",
            "CosmicExplorerPessimisticP1600143",
            "CosmicExplorerWidebandP1600143"};
ngLabels = {"ET-D", "CE (baseline)", "CE (pessimistic)", "CE (wideband)"};
ngColors = {Purple, Orange, Red, Green};

plotNG = Show[
  Table[
    With[
      {psd = DetectorPSD[ngModels[[i]], nFreq, deltaF, fLow],
       col = ngColors[[i]]},
      ListLogLogPlot[
        psdToASD[psd],
        Joined -> True,
        PlotStyle -> {col, Thick},
        PlotLegends -> None
      ]
    ],
    {i, Length[ngModels]}
  ],
  (* Overlay the aLIGO tabulated curve for reference *)
  ListLogLogPlot[asdGWINC, Joined -> True,
                 PlotStyle -> {Gray, Dashed, Thick}],
  PlotLabel -> "Next-Generation Ground-Based Detectors",
  AxesLabel -> {"f [Hz]", "ASD [1/\[Sqrt]Hz]"},
  PlotRange -> {{1, 4096}, All},
  GridLines -> Automatic, GridLinesStyle -> LightGray,
  PlotLegends -> Placed[
    LineLegend[
      Append[ngColors, Gray],
      Append[ngLabels, "aLIGO (reference)"],
      LegendMarkerSize -> 30
    ], {0.7, 0.3}],
  ImageSize -> 600
];

Export[FileNameJoin[{outputDir, "next_gen_detectors.pdf"}], plotNG];

(* ----------------------------------------------------------- *)
(* 4. LISA: effect of galactic confusion noise                 *)
(* ----------------------------------------------------------- *)

Print["-- LISA: galactic confusion noise --"];

lisaDeltaF = 1.0*^-4;  (* 0.1 mHz spacing for LISA band *)
lisaNFreq  = 40001;    (* up to 4 Hz *)
lisaFLow   = 1.0*^-4;

lisaFreqs = Table[(k - 1) * lisaDeltaF, {k, 2, lisaNFreq}];

lisaModels  = {"LISA", "LISAConfusion05yr", "LISAConfusion1yr",
               "LISAConfusion2yr", "LISAConfusion4yr"};
lisaLabels  = {"No confusion noise", "0.5 yr", "1 yr", "2 yr", "4 yr"};
lisaColors  = {Black, Blue, Cyan, Green, Red};

plotLISA = LogLogPlot[
  Evaluate[Table[Sqrt[($detectorPSDFunctions[m])[f]], {m, lisaModels}]],
  {f, 1.0*^-4, 0.5},
  PlotStyle -> MapThread[{#1, Thick} &, {lisaColors}],
  PlotLegends -> Placed[
    LineLegend[lisaColors, lisaLabels,
               LegendLabel -> "Galactic confusion noise\n(observation time)",
               LegendMarkerSize -> 30],
    {0.7, 0.3}],
  PlotLabel -> "LISA Sensitivity",
  AxesLabel -> {"f [Hz]", "ASD [1/\[Sqrt]Hz]"},
  PlotRange -> {{1*^-4, 0.5}, {1*^-20, 1*^-14}},
  GridLines -> Automatic, GridLinesStyle -> LightGray,
  ImageSize -> 600
];

Export[FileNameJoin[{outputDir, "LISA_confusion_noise.pdf"}], plotLISA];

(* ----------------------------------------------------------- *)
(* 5. Effective bandwidth for match calculations               *)
(* ----------------------------------------------------------- *)

Print["\n-- Effective bandwidth comparison (fLow=20 Hz, fHigh=1024 Hz) --"];

groundModels = {"aLIGOZeroDetHighPower", "AdvVirgo", "KAGRA",
                "EinsteinTelescopeP1600143", "CosmicExplorerP1600143"};

nFreqGround  = 2049;
deltaFGround = 0.5;
fLowGround   = 20.0;
fHighGround  = 1024.0;

Table[
  With[{psd = DetectorPSD[m, nFreqGround, deltaFGround, fLowGround],
        hTilde = ConstantArray[1.0, nFreqGround]},
    sigmasq = GWSigmaSq[hTilde, psd, deltaFGround, fLowGround, fHighGround];
    Print["  ", m, ":  sigma^2 = ", ScientificForm[sigmasq, 4]]
  ],
  {m, groundModels}
];

Print["\nPlots saved to: ", outputDir];
