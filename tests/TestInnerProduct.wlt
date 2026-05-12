(* TestInnerProduct.wlt *)
(* Unit tests for GWInnerProduct, GWSigmaSq, GWSigma, GWOverlap, *)
(* GetCutoffIndices, and CyclicTimeShift. *)
(* Run with: TestReport["path/to/TestInnerProduct.wlt"] *)

<< GWMMatch`

(* ------------------------------------------------------------------ *)
(* Shared fixtures *)
(* ------------------------------------------------------------------ *)

With[
  {
    (* Simple sinusoid at a single frequency *)
    sampleRate = 4096.0,
    duration   = 1.0,   (* seconds *)
    f0         = 100.0  (* Hz *)
  },

  (* Build a pure sinusoid time-domain signal *)
  $nTime  = Round[sampleRate * duration];
  $deltaT = 1.0 / sampleRate;
  $deltaF = 1.0 / duration;
  $fLow   = 10.0;
  $fHigh  = 1000.0;

  $h = Table[Sin[2.0 Pi f0 t], {t, 0, ($nTime - 1) $deltaT, $deltaT}];

  (* Flat PSD over the whole band *)
  $nFreq = Floor[$nTime / 2] + 1;
  $psd   = FlatPSD[$nFreq, $deltaF, $fLow];

  (* FFT of the sinusoid *)
  $hTilde = GWMMatch`Private`ForwardFFT[$h, $deltaT];
]

(* ------------------------------------------------------------------ *)
(* GetCutoffIndices *)
(* ------------------------------------------------------------------ *)

VerificationTest[
  (* kmin must be at least 2 (skip DC) *)
  First[GetCutoffIndices[10.0, 1000.0, 1.0, 4096]] >= 2,
  True,
  TestID -> "GetCutoffIndices_kmin_skips_DC"
]

VerificationTest[
  (* kmin = Floor[fLow/deltaF] + 1 *)
  First[GetCutoffIndices[20.0, 500.0, 1.0, 4096]],
  21,
  TestID -> "GetCutoffIndices_kmin_value"
]

VerificationTest[
  (* kmax = Floor[fHigh/deltaF] + 1 *)
  Last[GetCutoffIndices[20.0, 500.0, 1.0, 4096]],
  501,
  TestID -> "GetCutoffIndices_kmax_value"
]

VerificationTest[
  (* kmin < kmax always *)
  Apply[Less, GetCutoffIndices[10.0, 1000.0, 1.0, 4096]],
  True,
  TestID -> "GetCutoffIndices_kmin_lt_kmax"
]

(* ------------------------------------------------------------------ *)
(* GWSigmaSq / GWSigma *)
(* ------------------------------------------------------------------ *)

VerificationTest[
  (* sigmasq must be real and positive *)
  GWSigmaSq[$hTilde, $psd, $deltaF, $fLow, $fHigh] > 0,
  True,
  TestID -> "GWSigmaSq_positive"
]

VerificationTest[
  (* GWSigma = Sqrt[GWSigmaSq] *)
  Abs[GWSigma[$hTilde, $psd, $deltaF, $fLow, $fHigh]^2 -
      GWSigmaSq[$hTilde, $psd, $deltaF, $fLow, $fHigh]] < 1.0*^-10,
  True,
  TestID -> "GWSigma_equals_sqrt_sigmasq"
]

VerificationTest[
  (* Scaling h by a real constant c scales sigmasq by c^2 *)
  With[{c = 3.7},
    Abs[GWSigmaSq[c $hTilde, $psd, $deltaF, $fLow, $fHigh] /
        GWSigmaSq[$hTilde, $psd, $deltaF, $fLow, $fHigh] - c^2]
  ] < 1.0*^-10,
  True,
  TestID -> "GWSigmaSq_scales_with_amplitude_squared"
]

(* ------------------------------------------------------------------ *)
(* GWInnerProduct *)
(* ------------------------------------------------------------------ *)

VerificationTest[
  (* <h|h> must be real (imaginary part negligible) *)
  Abs[Im[GWInnerProduct[$hTilde, $hTilde, $psd, $deltaF, $fLow, $fHigh]]] < 1.0*^-10,
  True,
  TestID -> "GWInnerProduct_self_is_real"
]

VerificationTest[
  (* <h|h> equals GWSigmaSq[h] *)
  Abs[Re[GWInnerProduct[$hTilde, $hTilde, $psd, $deltaF, $fLow, $fHigh]] -
      GWSigmaSq[$hTilde, $psd, $deltaF, $fLow, $fHigh]] < 1.0*^-10,
  True,
  TestID -> "GWInnerProduct_self_equals_sigmasq"
]

VerificationTest[
  (* Conjugate symmetry: <h1|h2> = Conj[<h2|h1>] *)
  With[
    {h2Tilde = GWMMatch`Private`ForwardFFT[RotateLeft[$h, 50], $deltaT]},
    Abs[GWInnerProduct[$hTilde, h2Tilde, $psd, $deltaF, $fLow, $fHigh] -
        Conjugate[GWInnerProduct[h2Tilde, $hTilde, $psd, $deltaF, $fLow, $fHigh]]]
  ] < 1.0*^-10,
  True,
  TestID -> "GWInnerProduct_conjugate_symmetry"
]

VerificationTest[
  (* Linearity: <h | a*h> = a * <h|h> for real a *)
  With[{a = 2.5},
    Abs[GWInnerProduct[$hTilde, a $hTilde, $psd, $deltaF, $fLow, $fHigh] -
        a * GWInnerProduct[$hTilde, $hTilde, $psd, $deltaF, $fLow, $fHigh]]
  ] < 1.0*^-10,
  True,
  TestID -> "GWInnerProduct_linearity"
]

(* ------------------------------------------------------------------ *)
(* GWOverlap *)
(* ------------------------------------------------------------------ *)

VerificationTest[
  (* Self-overlap is exactly 1 *)
  Abs[GWOverlap[$hTilde, $hTilde, $psd, $deltaF, $fLow, $fHigh] - 1.0] < 1.0*^-10,
  True,
  TestID -> "GWOverlap_self_is_one"
]

VerificationTest[
  (* Overlap is in [-1, 1] *)
  With[{h2Tilde = GWMMatch`Private`ForwardFFT[Reverse[$h], $deltaT]},
    -1.0 <= GWOverlap[$hTilde, h2Tilde, $psd, $deltaF, $fLow, $fHigh] <= 1.0
  ],
  True,
  TestID -> "GWOverlap_in_range"
]

VerificationTest[
  (* Negating one waveform flips sign of the overlap *)
  Abs[GWOverlap[$hTilde, -$hTilde, $psd, $deltaF, $fLow, $fHigh] - (-1.0)] < 1.0*^-10,
  True,
  TestID -> "GWOverlap_negated_waveform"
]

(* ------------------------------------------------------------------ *)
(* CyclicTimeShift *)
(* ------------------------------------------------------------------ *)

VerificationTest[
  (* A zero time shift leaves the waveform unchanged *)
  Max[Abs[CyclicTimeShift[$hTilde, 0.0, $deltaF] - $hTilde]] < 1.0*^-12,
  True,
  TestID -> "CyclicTimeShift_zero_shift_identity"
]

VerificationTest[
  (* Time shift preserves the norm *)
  Abs[GWSigmaSq[CyclicTimeShift[$hTilde, 0.01, $deltaF],
                $psd, $deltaF, $fLow, $fHigh] -
      GWSigmaSq[$hTilde, $psd, $deltaF, $fLow, $fHigh]] < 1.0*^-10,
  True,
  TestID -> "CyclicTimeShift_preserves_norm"
]
