(* TestPSD.wlt *)
(* Unit tests for DetectorPSD, FlatPSD, LoadPSD, and ListDetectorPSDs. *)
(* Run with: TestReport["path/to/TestPSD.wlt"] *)

<< GWMMatch`

(* ------------------------------------------------------------------ *)
(* Shared parameters *)
(* ------------------------------------------------------------------ *)

$nFreq  = 4097;   (* 4096-sample FFT -> 4097 frequency bins *)
$deltaF = 1.0;    (* 1 Hz bin spacing *)
$fLow   = 20.0;

(* ------------------------------------------------------------------ *)
(* ListDetectorPSDs *)
(* ------------------------------------------------------------------ *)

VerificationTest[
  (* Returns a non-empty list of strings *)
  ListQ[ListDetectorPSDs[]] && Length[ListDetectorPSDs[]] > 0,
  True,
  TestID -> "ListDetectorPSDs_nonempty"
]

VerificationTest[
  (* All entries are strings *)
  AllTrue[ListDetectorPSDs[], StringQ],
  True,
  TestID -> "ListDetectorPSDs_all_strings"
]

VerificationTest[
  (* Key models are present *)
  SubsetQ[ListDetectorPSDs[], {"aLIGOZeroDetHighPower", "AdvVirgo", "KAGRA",
                               "EinsteinTelescopeP1600143", "LISA"}],
  True,
  TestID -> "ListDetectorPSDs_contains_key_models"
]

(* ------------------------------------------------------------------ *)
(* FlatPSD *)
(* ------------------------------------------------------------------ *)

VerificationTest[
  (* Length matches nFreq *)
  Length[FlatPSD[$nFreq, $deltaF, $fLow]],
  $nFreq,
  TestID -> "FlatPSD_correct_length"
]

VerificationTest[
  (* Bins below fLow are Infinity *)
  With[{psd = FlatPSD[$nFreq, $deltaF, $fLow]},
    AllTrue[psd[[1 ;; Floor[$fLow / $deltaF]]], # === Infinity &]
  ],
  True,
  TestID -> "FlatPSD_below_fLow_is_infinity"
]

VerificationTest[
  (* In-band values are 1.0 *)
  With[{psd = FlatPSD[$nFreq, $deltaF, $fLow]},
    AllTrue[psd[[Floor[$fLow / $deltaF] + 2 ;; -1]], # == 1.0 &]
  ],
  True,
  TestID -> "FlatPSD_inband_is_unity"
]

(* ------------------------------------------------------------------ *)
(* DetectorPSD: general contracts for all analytical models *)
(* ------------------------------------------------------------------ *)

(* Analytical models that do not require data files *)
$analyticalModels = {
  "aLIGOZeroDetHighPower", "aLIGOZeroDetLowPower", "aLIGONSNSOpt",
  "aLIGOBHBH20Deg", "aLIGOHighFrequency",
  "AdvVirgo", "KAGRA", "Virgo", "iLIGOSRD",
  "LISA", "LISAConfusion05yr", "LISAConfusion1yr",
  "LISAConfusion2yr", "LISAConfusion4yr"
};

Do[
  With[{name = model, psd = DetectorPSD[model, $nFreq, $deltaF, $fLow]},
    VerificationTest[
      Length[psd] == $nFreq,
      True,
      TestID -> "DetectorPSD_length_" <> name
    ];
    VerificationTest[
      (* Bins at DC and below fLow are Infinity *)
      psd[[1]] === Infinity,
      True,
      TestID -> "DetectorPSD_DC_is_infinity_" <> name
    ];
    VerificationTest[
      (* In-band values are finite and positive *)
      With[{inBand = Select[psd, # =!= Infinity &]},
        AllTrue[inBand, # > 0 &]
      ],
      True,
      TestID -> "DetectorPSD_inband_positive_" <> name
    ]
  ],
  {model, $analyticalModels}
]

(* ------------------------------------------------------------------ *)
(* DetectorPSD: tabulated models (require psd_data/ files) *)
(* ------------------------------------------------------------------ *)

$tabulatedModels = {
  "aLIGOZeroDetHighPowerGWINC",
  "EinsteinTelescopeP1600143",
  "CosmicExplorerP1600143",
  "CosmicExplorerPessimisticP1600143",
  "CosmicExplorerWidebandP1600143"
};

Do[
  With[{name = model, psd = DetectorPSD[model, $nFreq, $deltaF, $fLow]},
    VerificationTest[
      psd =!= $Failed,
      True,
      TestID -> "DetectorPSD_tabulated_loaded_" <> name
    ];
    VerificationTest[
      Length[psd] == $nFreq,
      True,
      TestID -> "DetectorPSD_tabulated_length_" <> name
    ];
    VerificationTest[
      With[{inBand = Select[psd, # =!= Infinity &]},
        AllTrue[inBand, # > 0 &]
      ],
      True,
      TestID -> "DetectorPSD_tabulated_inband_positive_" <> name
    ]
  ],
  {model, $tabulatedModels}
]

(* ------------------------------------------------------------------ *)
(* DetectorPSD: unknown model returns $Failed *)
(* ------------------------------------------------------------------ *)

VerificationTest[
  DetectorPSD["NonExistentDetector", $nFreq, $deltaF, $fLow],
  $Failed,
  TestID -> "DetectorPSD_unknown_model_returns_failed"
]

(* ------------------------------------------------------------------ *)
(* LoadPSD *)
(* ------------------------------------------------------------------ *)

VerificationTest[
  (* LoadPSD from the aLIGO GWINC file *)
  With[
    {filename = FileNameJoin[{
       DirectoryName[FindFile["GWMMatch`"]],
       "psd_data", "LIGO-T0900288-v3-ZERO_DET_high_P.txt"}],
     psd = LoadPSD[
       FileNameJoin[{DirectoryName[FindFile["GWMMatch`"]],
                     "psd_data", "LIGO-T0900288-v3-ZERO_DET_high_P.txt"}],
       $deltaF, $fLow, $nFreq]},
    ListQ[psd] && Length[psd] == $nFreq
  ],
  True,
  TestID -> "LoadPSD_correct_length"
]

VerificationTest[
  (* LoadPSD: DC bin is Infinity *)
  With[
    {psd = LoadPSD[
       FileNameJoin[{DirectoryName[FindFile["GWMMatch`"]],
                     "psd_data", "LIGO-T0900288-v3-ZERO_DET_high_P.txt"}],
       $deltaF, $fLow, $nFreq]},
    psd[[1]] === Infinity
  ],
  True,
  TestID -> "LoadPSD_DC_is_infinity"
]

VerificationTest[
  (* LoadPSD: in-band values are finite and positive *)
  With[
    {psd = LoadPSD[
       FileNameJoin[{DirectoryName[FindFile["GWMMatch`"]],
                     "psd_data", "LIGO-T0900288-v3-ZERO_DET_high_P.txt"}],
       $deltaF, $fLow, $nFreq]},
    AllTrue[Select[psd, # =!= Infinity &], # > 0 &]
  ],
  True,
  TestID -> "LoadPSD_inband_positive"
]
