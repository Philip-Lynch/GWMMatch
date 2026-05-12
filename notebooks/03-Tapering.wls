(* =========================================================== *)
(*  03-Tapering.wls                                            *)
(*  GWMMatch — Planck tapering: time and frequency domain     *)
(*                                                             *)
(*  Run:  wolframscript -file 03-Tapering.wls                 *)
(*  Or:   Get["path/to/notebooks/03-Tapering.wls"]           *)
(*                                                             *)
(*  Plots are saved to 03-Tapering-output/ as PDF files.      *)
(* =========================================================== *)

<< GWMMatch`

outputDir = FileNameJoin[{DirectoryName[$InputFileName], "03-Tapering-output"}];
If[!DirectoryQ[outputDir], CreateDirectory[outputDir]];

(* ----------------------------------------------------------- *)
(* 1. Build a test signal                                      *)
(* ----------------------------------------------------------- *)

sampleRate = 4096.0;
duration   = 1.0;
nTime      = Round[sampleRate * duration];
deltaT     = 1.0 / sampleRate;
deltaF     = 1.0 / duration;
fLow       = 20.0;
fHigh      = 1000.0;

(* A pure sinusoid with a sharp onset — pathological spectral leakage *)
f0 = 150.0;  (* Hz *)
tVals  = Table[n * deltaT, {n, 0, nTime - 1}];
hSharp = Table[Sin[2.0 Pi f0 t], {t, tVals}];

(* Apply tapers with three different widths *)
hTaperAuto  = GWTaper[hSharp];                             (* auto-detect width *)
hTaperFixed = GWTaper[hSharp, "TAPER_STARTEND", 200];     (* 200 samples each end *)
hTaperStart = GWTaper[hSharp, "TAPER_START"];              (* start only *)

Print["Signal length: ", nTime, " samples"];

(* ----------------------------------------------------------- *)
(* 2. Time-domain view                                         *)
(* ----------------------------------------------------------- *)

nShow = 400;   (* samples to display at each edge *)

plotTimeStart = ListPlot[
  {
    Transpose[{tVals[[1 ;; nShow]], hSharp[[1 ;; nShow]]}],
    Transpose[{tVals[[1 ;; nShow]], hTaperAuto[[1 ;; nShow]]}],
    Transpose[{tVals[[1 ;; nShow]], hTaperFixed[[1 ;; nShow]]}]
  },
  Joined -> True,
  PlotStyle -> {{Blue, Thin}, {Red, Thick}, {Green, Thick, Dashed}},
  PlotLegends -> {"Untapered", "Planck taper (auto)", "Planck taper (200 samples)"},
  PlotLabel -> "Waveform start — effect of Planck tapering",
  AxesLabel -> {"t [s]", "h(t)"},
  PlotRange -> All, ImageSize -> 600
];

plotTimeEnd = ListPlot[
  {
    Transpose[{tVals[[-nShow ;; -1]], hSharp[[-nShow ;; -1]]}],
    Transpose[{tVals[[-nShow ;; -1]], hTaperAuto[[-nShow ;; -1]]}],
    Transpose[{tVals[[-nShow ;; -1]], hTaperFixed[[-nShow ;; -1]]}]
  },
  Joined -> True,
  PlotStyle -> {{Blue, Thin}, {Red, Thick}, {Green, Thick, Dashed}},
  PlotLegends -> {"Untapered", "Planck taper (auto)", "Planck taper (200 samples)"},
  PlotLabel -> "Waveform end — effect of Planck tapering",
  AxesLabel -> {"t [s]", "h(t)"},
  PlotRange -> All, ImageSize -> 600
];

Export[FileNameJoin[{outputDir, "taper_time_start.pdf"}], plotTimeStart];
Export[FileNameJoin[{outputDir, "taper_time_end.pdf"}], plotTimeEnd];
Print["Time-domain plots saved."];

(* ----------------------------------------------------------- *)
(* 3. Frequency-domain: spectral leakage                      *)
(* ----------------------------------------------------------- *)

(* Compute one-sided ASDs *)
asdOf[h_] := With[
  {nF = Floor[nTime / 2] + 1,
   fft = GWMMatch`Private`ForwardFFT[h, deltaT]},
  Transpose[{
    Table[(k - 1) * deltaF, {k, 1, nF}],
    Abs[fft] * 2.0 * deltaT  (* one-sided ASD, unnormalised *)
  }]
];

asdSharp = asdOf[hSharp];
asdAuto  = asdOf[hTaperAuto];
asdFixed = asdOf[hTaperFixed];

(* Zoom into the spectral leakage region around f0 *)
fBand = {100, 200};  (* Hz *)

plotASD = ListLogLogPlot[
  {asdSharp, asdAuto, asdFixed},
  Joined -> True,
  PlotStyle -> {{Blue, Thin}, {Red, Thick}, {Green, Thick, Dashed}},
  PlotLegends -> {"Untapered", "Planck taper (auto)", "Planck taper (200 samples)"},
  PlotLabel -> "Spectral leakage around " <> ToString[f0] <> " Hz",
  AxesLabel -> {"f [Hz]", "|h(f)|"},
  PlotRange -> {fBand, All},
  GridLines -> Automatic, GridLinesStyle -> LightGray,
  ImageSize -> 600
];

plotASDFull = ListLogLogPlot[
  {asdSharp, asdAuto, asdFixed},
  Joined -> True,
  PlotStyle -> {{Blue, Thin}, {Red, Thick}, {Green, Thick, Dashed}},
  PlotLegends -> {"Untapered", "Planck taper (auto)", "Planck taper (200 samples)"},
  PlotLabel -> "Full spectrum — spectral leakage suppression",
  AxesLabel -> {"f [Hz]", "|h(f)|"},
  PlotRange -> {{10, 2000}, All},
  GridLines -> Automatic, GridLinesStyle -> LightGray,
  ImageSize -> 600
];

Export[FileNameJoin[{outputDir, "taper_spectral_zoom.pdf"}], plotASD];
Export[FileNameJoin[{outputDir, "taper_spectral_full.pdf"}], plotASDFull];
Print["Frequency-domain plots saved."];

(* ----------------------------------------------------------- *)
(* 4. Effect of tapering on match                             *)
(* ----------------------------------------------------------- *)

nFreq   = Floor[nTime / 2] + 1;
psdFlat = FlatPSD[nFreq, deltaF, fLow];

(* A slightly different sinusoid at 151 Hz *)
hRef = Table[Sin[2.0 Pi 151.0 t], {t, tVals}];
hRefTapered = GWTaper[hRef];

computeMatch[ha_, hb_] :=
  First[GWOptimizedMatch[ha, hb, deltaT, psdFlat, deltaF, fLow, fHigh]];

Print["\n-- Match between 150 Hz and 151 Hz sinusoids --"];
Print["  Untapered vs untapered:  ", N[computeMatch[hSharp, hRef], 6]];
Print["  Both tapered (auto):     ", N[computeMatch[hTaperAuto, hRefTapered], 6]];
Print["  Both tapered (200 samp): ", N[computeMatch[hTaperFixed, GWTaper[hRef, "TAPER_STARTEND", 200]], 6]];

(* ----------------------------------------------------------- *)
(* 5. Taper width parameter study                             *)
(* ----------------------------------------------------------- *)

Print["\n-- Match vs taper width (self-match of tapered vs untapered) --"];

widths = {20, 50, 100, 200, 400};

Table[
  With[{ht = GWTaper[hSharp, "TAPER_STARTEND", w]},
    m = N[computeMatch[hSharp, ht], 6];
    Print["  width = ", w, " samples: match = ", m]
  ],
  {w, widths}
];

Print["\nAll outputs saved to: ", outputDir];
