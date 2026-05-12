(* TestMatch.wlt *)
(* Unit tests for GWMatch, GWOptimizedMatch, GWMismatch, GWTaper. *)
(* Run with: TestReport["path/to/TestMatch.wlt"] *)

<< GWMMatch`

(* ------------------------------------------------------------------ *)
(* Shared fixtures *)
(* ------------------------------------------------------------------ *)

With[
  {
    sampleRate = 2048.0,
    duration   = 2.0,
    f0         = 60.0
  },

  $nTime  = Round[sampleRate * duration];
  $deltaT = 1.0 / sampleRate;
  $deltaF = 1.0 / duration;
  $fLow   = 20.0;
  $fHigh  = 512.0;
  $nFreq  = Floor[$nTime / 2] + 1;
  $psd    = FlatPSD[$nFreq, $deltaF, $fLow];

  (* A smooth, band-limited chirp-like signal *)
  $h1 = Table[
    Sin[2.0 Pi (f0 + 5.0 t) t] * Exp[-(t - duration/2)^2 / 0.1],
    {t, 0, ($nTime - 1) $deltaT, $deltaT}
  ];

  (* h2: same signal shifted by 50 samples *)
  $shiftSamples = 50;
  $h2 = RotateLeft[$h1, $shiftSamples];
]

(* ------------------------------------------------------------------ *)
(* GWMatch: identical waveforms *)
(* ------------------------------------------------------------------ *)

VerificationTest[
  (* match(h, h) = 1 *)
  Abs[First[GWMatch[$h1, $h1, $deltaT, $psd, $deltaF, $fLow, $fHigh]] - 1.0] < 1.0*^-6,
  True,
  TestID -> "GWMatch_identical_waveforms"
]

VerificationTest[
  (* match(h, -h) = 1 — phase maximization absorbs the sign flip *)
  Abs[First[GWMatch[$h1, -$h1, $deltaT, $psd, $deltaF, $fLow, $fHigh]] - 1.0] < 1.0*^-6,
  True,
  TestID -> "GWMatch_negated_waveform"
]

VerificationTest[
  (* match(h, A*h) = 1 for any positive amplitude A *)
  Abs[First[GWMatch[$h1, 7.3 $h1, $deltaT, $psd, $deltaF, $fLow, $fHigh]] - 1.0] < 1.0*^-6,
  True,
  TestID -> "GWMatch_amplitude_invariance"
]

VerificationTest[
  (* match(h1, shifted h1) = 1 — time maximization recovers the shift *)
  Abs[First[GWMatch[$h1, $h2, $deltaT, $psd, $deltaF, $fLow, $fHigh]] - 1.0] < 1.0*^-4,
  True,
  TestID -> "GWMatch_time_shifted_waveform"
]

VerificationTest[
  (* match is in [0, 1] *)
  0.0 <= First[GWMatch[$h1, $h2, $deltaT, $psd, $deltaF, $fLow, $fHigh]] <= 1.0,
  True,
  TestID -> "GWMatch_range"
]

VerificationTest[
  (* GWMatch returns a list of length 2: {match, timeIndex} *)
  Length[GWMatch[$h1, $h2, $deltaT, $psd, $deltaF, $fLow, $fHigh]],
  2,
  TestID -> "GWMatch_return_shape"
]

VerificationTest[
  (* Length mismatch returns $Failed *)
  GWMatch[$h1, Most[$h1], $deltaT, $psd, $deltaF, $fLow, $fHigh],
  $Failed,
  TestID -> "GWMatch_length_mismatch"
]

(* ------------------------------------------------------------------ *)
(* GWOptimizedMatch *)
(* ------------------------------------------------------------------ *)

VerificationTest[
  (* optimized match(h, h) = 1 *)
  Abs[First[GWOptimizedMatch[$h1, $h1, $deltaT, $psd, $deltaF, $fLow, $fHigh]] - 1.0] < 1.0*^-6,
  True,
  TestID -> "GWOptimizedMatch_identical_waveforms"
]

VerificationTest[
  (* optimized match(h, -h) = 1 *)
  Abs[First[GWOptimizedMatch[$h1, -$h1, $deltaT, $psd, $deltaF, $fLow, $fHigh]] - 1.0] < 1.0*^-6,
  True,
  TestID -> "GWOptimizedMatch_negated_waveform"
]

VerificationTest[
  (* optimized match >= discrete match (it can only be better or equal) *)
  First[GWOptimizedMatch[$h1, $h2, $deltaT, $psd, $deltaF, $fLow, $fHigh]] >=
    First[GWMatch[$h1, $h2, $deltaT, $psd, $deltaF, $fLow, $fHigh]] - 1.0*^-6,
  True,
  TestID -> "GWOptimizedMatch_ge_discrete_match"
]

VerificationTest[
  (* optimized match(h1, shifted h1) is close to 1 *)
  Abs[First[GWOptimizedMatch[$h1, $h2, $deltaT, $psd, $deltaF, $fLow, $fHigh]] - 1.0] < 1.0*^-5,
  True,
  TestID -> "GWOptimizedMatch_time_shifted_waveform"
]

VerificationTest[
  (* GWOptimizedMatch returns a list of length 2 *)
  Length[GWOptimizedMatch[$h1, $h2, $deltaT, $psd, $deltaF, $fLow, $fHigh]],
  2,
  TestID -> "GWOptimizedMatch_return_shape"
]

(* ------------------------------------------------------------------ *)
(* GWMismatch *)
(* ------------------------------------------------------------------ *)

VerificationTest[
  (* mismatch(h, h) = 0 *)
  Abs[First[GWMismatch[$h1, $h1, $deltaT, $psd, $deltaF, $fLow, $fHigh]]] < 1.0*^-6,
  True,
  TestID -> "GWMismatch_identical_waveforms"
]

VerificationTest[
  (* mismatch = 1 - match *)
  With[
    {mm = First[GWMismatch[$h1, $h2, $deltaT, $psd, $deltaF, $fLow, $fHigh]],
     m  = First[GWOptimizedMatch[$h1, $h2, $deltaT, $psd, $deltaF, $fLow, $fHigh]]},
    Abs[mm - (1.0 - m)]
  ] < 1.0*^-10,
  True,
  TestID -> "GWMismatch_equals_one_minus_match"
]

VerificationTest[
  (* mismatch is in [0, 1] *)
  0.0 <= First[GWMismatch[$h1, $h2, $deltaT, $psd, $deltaF, $fLow, $fHigh]] <= 1.0,
  True,
  TestID -> "GWMismatch_range"
]

(* ------------------------------------------------------------------ *)
(* GWTaper *)
(* ------------------------------------------------------------------ *)

With[
  {
    nTaper   = 4096,
    hSimple  = Join[ConstantArray[0.0, 100],
                    Table[Sin[2 Pi 50.0 t / 2048.0], {t, 0, 3795}],
                    ConstantArray[0.0, 100]]
  },

  VerificationTest[
    (* Output has the same length as input *)
    Length[GWTaper[hSimple]] == Length[hSimple],
    True,
    TestID -> "GWTaper_preserves_length"
  ];

  VerificationTest[
    (* First nonzero sample after a start taper should be zero (Planck sets endpoints = 0) *)
    With[{tapered = GWTaper[hSimple, "TAPER_START"]},
      (* Find the first nonzero input sample and confirm it is zeroed *)
      Module[{start = First[FirstPosition[hSimple, x_ /; x != 0.0]]},
        tapered[[start]] == 0.0
      ]
    ],
    True,
    TestID -> "GWTaper_start_endpoint_is_zero"
  ];

  VerificationTest[
    (* Last nonzero sample after an end taper should be zero *)
    With[{tapered = GWTaper[hSimple, "TAPER_END"]},
      Module[{endd = First[First[Position[hSimple, x_ /; x != 0.0, 1, 1, Heads -> False], 1] /.
                First -> Last]},
        (* Use Last position with nonzero value *)
        Module[{lastIdx = Last[Flatten[Position[hSimple, x_ /; x != 0.0]]]},
          tapered[[lastIdx]] == 0.0
        ]
      ]
    ],
    True,
    TestID -> "GWTaper_end_endpoint_is_zero"
  ];

  VerificationTest[
    (* Fixed taper width is respected: nEnd samples from end are modified *)
    With[
      {hFlat = Join[ConstantArray[0.0, 50],
                    Table[1.0, {3900}],
                    ConstantArray[0.0, 50]],
       nFixed = 100},
      With[
        {tapered = GWTaper[hFlat, "TAPER_STARTEND", nFixed]},
        (* Interior of the flat section should be untouched at 1.0 *)
        tapered[[200]] == 1.0
      ]
    ],
    True,
    TestID -> "GWTaper_interior_untouched"
  ]
]
