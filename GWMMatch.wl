(* ::Package:: *)
(* GWMMatch.wl *)
(* Gravitational Wave Match/Mismatch Package *)
(* Port of pyCBC match/optimized_match to Mathematica *)
(* Reference: pycbc.filter.matchedfilter (Alex Nitz et al.) *)

BeginPackage["GWMMatch`"];

(* ===================== Public Interface ===================== *)

GWInnerProduct::usage = 
"GWInnerProduct[h1Tilde, h2Tilde, psd, deltaF, fLow, fHigh] computes the \
noise-weighted inner product <h1|h2> = 4 Re[Sum[Conj[h1]*h2/Sn * deltaF]] \
over the frequency band [fLow, fHigh]. h1Tilde, h2Tilde, and psd are lists \
of complex/real values sampled at frequencies k*deltaF. Returns a complex number; \
take Re[] for the real inner product.";

GWSigmaSq::usage = 
"GWSigmaSq[hTilde, psd, deltaF, fLow, fHigh] computes the sigmasq = <h|h> \
= 4 Sum[|h(f)|^2 / Sn(f)] * deltaF over [fLow, fHigh].";

GWSigma::usage = 
"GWSigma[hTilde, psd, deltaF, fLow, fHigh] returns Sqrt[GWSigmaSq[...]].";

GWOverlap::usage = 
"GWOverlap[h1Tilde, h2Tilde, psd, deltaF, fLow, fHigh] computes the \
normalized overlap (no maximization over time/phase): \
Re[<h1|h2>] / (sigma1 * sigma2). Returns a real number in [-1,1].";

GWMatch::usage = 
"GWMatch[h1, h2, deltaT, psd, deltaF, fLow, fHigh, opts] computes the match \
between two time-domain waveforms h1 and h2 (lists of real values sampled at \
interval deltaT), maximized over time and phase shifts. The PSD is given as a \
list sampled at spacing deltaF. Returns {match, timeIndex} where match is in \
[0,1]. Options: \"SubsampleInterpolation\" -> True/False (default True).";

GWOptimizedMatch::usage = 
"GWOptimizedMatch[h1, h2, deltaT, psd, deltaF, fLow, fHigh] computes the \
match between two time-domain waveforms, using Brent's method to find the \
optimal sub-sample time shift. This is more accurate than GWMatch and matches \
pyCBC's optimized_match function. Returns {match, timeShiftIndex}.";

GWMismatch::usage = 
"GWMismatch[h1, h2, deltaT, psd, deltaF, fLow, fHigh] returns \
1 - GWOptimizedMatch[...]. Convenience wrapper.";

LoadPSD::usage = 
"LoadPSD[filename, deltaF, fLow, nFreq] loads a two-column (frequency, PSD) \
text file and interpolates it onto a regular grid of nFreq points with spacing \
deltaF, starting at f=0. Values below fLow are set to Infinity (excluded from \
inner products). Returns a list of length nFreq.";

FlatPSD::usage = 
"FlatPSD[nFreq, deltaF, fLow] returns a flat (unity) PSD of length nFreq \
with values below fLow set to Infinity.";

CyclicTimeShift::usage = 
"CyclicTimeShift[hTilde, deltaT, deltaF, nFreq] applies a cyclic time shift \
of deltaT seconds to a frequency-domain waveform hTilde. Equivalent to \
multiplying by Exp[2 Pi I f deltaT].";

GetCutoffIndices::usage = 
"GetCutoffIndices[fLow, fHigh, deltaF, nTime] returns {kmin, kmax} frequency \
bin indices (1-based) corresponding to the frequency band [fLow, fHigh]. \
Matches pyCBC's get_cutoff_indices (converted to 1-based indexing).";

DetectorPSD::usage = 
"DetectorPSD[name, nFreq, deltaF, fLow] returns a one-sided PSD array of \
length nFreq sampled at spacing deltaF for the named detector. Bins below \
fLow are set to Infinity. Available names include analytical models \
(\"aLIGOZeroDetHighPower\", etc.) and tabulated models from LIGO-P1600143 \
(\"aLIGOZeroDetHighPowerGWINC\", \"EinsteinTelescopeP1600143\", \
\"CosmicExplorerP1600143\", \"CosmicExplorerPessimisticP1600143\", \
\"CosmicExplorerWidebandP1600143\") and space-based models \
(\"LISA\", \"LISAConfusion4yr\", etc.). \
Use ListDetectorPSDs[] for the full list.";

GWTaper::usage = 
"GWTaper[h, taperMethod, taperSamples] applies a Planck taper to the \
time-domain waveform h, matching LALSuite's XLALSimInspiralREAL8WaveTaper. \
taperMethod can be \"TAPER_START\", \"TAPER_END\", or \"TAPER_STARTEND\" \
(default). taperSamples controls the taper width: Automatic (default) \
auto-detects from waveform extrema; a single integer n uses that width at \
both ends; {nStart, nEnd} sets each end independently, where either element \
can be Automatic. Examples: GWTaper[h] auto-detects both ends; \
GWTaper[h, \"TAPER_STARTEND\", {Automatic, 200}] auto-detects the start \
but uses 200 samples for the end taper.";

ListDetectorPSDs::usage = 
"ListDetectorPSDs[] returns a list of all available detector PSD model names.";

Begin["`Private`"];

(* ===================== Utility Functions ===================== *)

(* Convert frequency cutoffs to 1-based indices into a frequency series.
   pyCBC uses 0-based; we add 1 for Mathematica. 
   nTime = number of time-domain samples (even). *)
GetCutoffIndices[fLow_, fHigh_, deltaF_, nTime_] := Module[
  {kmin, kmax},
  (* kmin: 1-based index *)
  If[fLow === None || fLow === 0,
    kmin = 2, (* skip DC, pyCBC uses kmin=1 in 0-based = index 2 in 1-based *)
    kmin = Floor[fLow / deltaF] + 1;
    If[kmin < 1, Message[GetCutoffIndices::negfreq]; Abort[]]
  ];
  (* kmax: 1-based index (exclusive upper bound like Python) *)
  If[fHigh === None || fHigh === 0,
    kmax = Floor[(nTime + 1) / 2] + 1, (* Nyquist, 1-based *)
    kmax = Floor[fHigh / deltaF] + 1;
    If[kmax > Floor[(nTime + 1) / 2] + 1,
      kmax = Floor[(nTime + 1) / 2] + 1
    ]
  ];
  If[kmax <= kmin, 
    Message[GetCutoffIndices::badrange]; Abort[]
  ];
  {kmin, kmax}
];

(* Frequency array for a frequency series of length nFreq with spacing deltaF *)
FrequencyArray[nFreq_, deltaF_] := Table[(k - 1) * deltaF, {k, 1, nFreq}];

(* ===================== Planck Taper ===================== *)

(* Planck taper sigmoid: sigma(x, n) = 1/(Exp[z] + 1)
   where z = (n-1)/x + (n-1)/(x - (n-1)).
   Faithful port of XLALSimInspiralREAL8WaveTaper from LALSuite.
   Reference: McKechan, Robinson, Sathyaprakash, CQG 27 (2010) 084020,
   Eq. 3.35 of gr-qc/0001023. *)
planckSigma[x_, n_] := Module[{z},
  z = (n - 1.0) / x + (n - 1.0) / (x - (n - 1.0));
  1.0 / (Exp[z] + 1.0)
];

(* Find the index of the 2nd valid local extremum of |signal| from the start.
   "Valid" means more than ringingExtent samples from the boundary.
   Returns the distance n from start to the 2nd peak. *)
findTaperN[data_List, start_Integer, mid_Integer, direction_String] := Module[
  {n, flag = 0, i, absData, ringingExtent = 19, len = Length[data]},
  absData = Abs[data];
  
  If[direction === "forward",
    n = mid - start;
    i = start + 1;
    While[i < mid,
      (* Check for local maximum of |data| *)
      If[absData[[i]] >= absData[[i - 1]] && absData[[i]] >= absData[[i + 1]],
        (* Skip flat peaks *)
        While[i + 1 <= mid && absData[[i]] == absData[[i + 1]], i++];
        n = i - start;
        If[i - start > ringingExtent, flag++];
        If[flag >= 2, Break[]];
      ];
      i++;
    ];
    n
    ,
    (* direction === "backward" *)
    n = start - mid;  (* start is the 'end' index here *)
    i = start - 1;
    While[i > mid,
      If[absData[[i]] >= absData[[i + 1]] && absData[[i]] >= absData[[i - 1]],
        While[i - 1 >= mid && absData[[i]] == absData[[i - 1]], i--];
        n = start - i;
        If[start - i > ringingExtent, flag++];
        If[flag >= 2, Break[]];
      ];
      i--;
    ];
    n
  ]
];

(* Apply Planck taper to a waveform. Matches XLALSimInspiralREAL8WaveTaper.
   Optional third argument: taper width in samples.
   Automatic = auto-detect from extrema (LALSuite behavior).
   Integer n = fixed width n at both ends.
   {nStart, nEnd} = independent widths; either can be Automatic. *)
GWTaper[hIn_List, taperMethod_String: "TAPER_STARTEND", taperSamples_: Automatic] := Module[
  {h = N[hIn], len, startIdx, endIdx, mid, n, x, nStart, nEnd,
   startSpec, endSpec},
  len = Length[h];
  
  (* Find first and last nonzero samples (1-based) *)
  startIdx = 1;
  While[startIdx <= len && h[[startIdx]] == 0.0, startIdx++];
  endIdx = len;
  While[endIdx >= 1 && h[[endIdx]] == 0.0, endIdx--];
  
  If[endIdx - startIdx <= 1, Return[h]];
  mid = Floor[(startIdx + endIdx) / 2];
  
  (* Parse taperSamples into per-end specs *)
  If[ListQ[taperSamples],
    startSpec = taperSamples[[1]];
    endSpec = taperSamples[[2]];
    ,
    startSpec = taperSamples;
    endSpec = taperSamples;
  ];
  
  (* Determine start taper width *)
  If[taperMethod === "TAPER_START" || taperMethod === "TAPER_STARTEND",
    nStart = If[startSpec === Automatic,
      findTaperN[h, startIdx, mid, "forward"],
      Min[startSpec, mid - startIdx]
    ];
    ,
    nStart = 0;
  ];
  
  (* Determine end taper width *)
  If[taperMethod === "TAPER_END" || taperMethod === "TAPER_STARTEND",
    nEnd = If[endSpec === Automatic,
      findTaperN[h, endIdx, mid, "backward"],
      Min[endSpec, endIdx - mid]
    ];
    ,
    nEnd = 0;
  ];
  
  (* Taper the start *)
  If[nStart > 1,
    h[[startIdx]] = 0.0;
    Do[
      x = i - startIdx;
      h[[i]] *= planckSigma[N[x], N[nStart]];
      , {i, startIdx + 1, startIdx + nStart - 2}
    ];
  ];
  
  (* Taper the end *)
  If[nEnd > 1,
    h[[endIdx]] = 0.0;
    Do[
      x = endIdx - i;
      h[[i]] *= planckSigma[N[x], N[nEnd]];
      , {i, endIdx - 1, endIdx - nEnd + 2, -1}
    ];
  ];
  
  h
];

(* ===================== FFT Conventions ===================== *)

(* pyCBC FFT convention:
   Forward:  hTilde[k] = deltaT * Sum[h[j] * Exp[-2 Pi I j k / N], {j,0,N-1}]
   Inverse (c2c): h[j] = Sum[hTilde[k] * Exp[2 Pi I j k / N], {k,0,N-1}]  (unnormalized)
   
   Mathematica Fourier with FourierParameters -> {1, -1}:
   F[k] = Sum[h[j] * Exp[-2 Pi I (j-1)(k-1) / N], {j,1,N}]  (unnormalized DFT)
   
   Mathematica InverseFourier with FourierParameters -> {-1, -1}:
   f[r] = Sum[F[s] * Exp[+2 Pi I (s-1)(r-1) / N], {s,1,N}]  (unnormalized IDFT)
   
   pyCBC matched_filter_core does a complex-to-complex IFFT on a full 
   length-N array where only positive frequency bins [kmin,kmax) are filled.
   The negative frequency bins remain zero. This produces a COMPLEX SNR 
   time series whose absolute value maximizes over phase. *)

(* Forward FFT: time domain (length N, spacing deltaT) -> frequency domain (length Floor[N/2]+1) *)
ForwardFFT[h_List, deltaT_] := Module[
  {n = Length[h], fullFFT, nFreq},
  nFreq = Floor[n / 2] + 1;
  fullFFT = Fourier[h, FourierParameters -> {1, -1}];
  deltaT * fullFFT[[1 ;; nFreq]]
];

(* c2c Inverse FFT on a full length-N complex array (unnormalized IDFT).
   This matches pyCBC's ifft(Array, Array) which uses FFTW backward with no 1/N. *)
InverseFFTc2c[fullSpectrum_List] := 
  InverseFourier[fullSpectrum, FourierParameters -> {-1, -1}];

(* ===================== Brent's Method ===================== *)

(* Brent's method for bounded 1D minimization.
   Faithful port of scipy.optimize._optimize._minimize_scalar_bounded.
   Variable names match scipy: xf=best x, nfc=2nd best, fulc=3rd best.
   Returns {fmin, xmin}. *)
BrentMinimize[f_, {x1In_?NumericQ, x2In_?NumericQ}, xatol_: 1.0*^-12, maxiter_: 500] := 
 Module[
  {a = N[x1In], b = N[x2In], 
   goldenMean = 0.5*(3.0 - Sqrt[5.0]),
   sqrtEps = Sqrt[$MachineEpsilon],
   fulc, nfc, xf, rat, e, fx, ffulc, fnfc,
   xm, tol1, tol2, golden, r, q, p, si, x, fu,
   num = 1},
  
  (* Initialize *)
  fulc = a + goldenMean*(b - a);
  nfc = xf = fulc;
  rat = e = 0.0;
  x = xf;
  fx = f[x];
  ffulc = fnfc = fx;
  xm = 0.5*(a + b);
  tol1 = sqrtEps*Abs[xf] + xatol/3.0;
  tol2 = 2.0*tol1;
  
  While[Abs[xf - xm] > (tol2 - 0.5*(b - a)) && num < maxiter,
    golden = True;
    
    If[Abs[e] > tol1,
      golden = False;
      (* Parabolic interpolation *)
      r = (xf - nfc)*(fx - ffulc);
      q = (xf - fulc)*(fx - fnfc);
      p = (xf - fulc)*q - (xf - nfc)*r;
      q = 2.0*(q - r);
      If[q > 0.0, p = -p];
      q = Abs[q];
      (* Key: save old e, then set e = previous rat *)
      r = e;
      e = rat;
      
      If[Abs[p] < Abs[0.5*q*r] && p > q*(a - xf) && p < q*(b - xf),
        (* Accept parabolic step *)
        rat = p/q;
        x = xf + rat;
        If[(x - a) < tol2 || (b - x) < tol2,
          si = Sign[xm - xf] + If[xm - xf == 0, 1, 0];
          rat = tol1*si;
        ],
        (* Reject parabolic, use golden section *)
        golden = True
      ]
    ];
    
    If[golden,
      If[xf >= xm,
        e = a - xf,
        e = b - xf
      ];
      rat = goldenMean*e;
    ];
    
    (* Evaluate function at new point *)
    si = Sign[rat] + If[rat == 0, 1, 0];
    x = xf + Max[Abs[rat], tol1]*si;
    fu = f[x];
    num++;
    
    (* Update bracket and best points *)
    If[fu <= fx,
      If[x >= xf, a = xf, b = xf];
      fulc = nfc; ffulc = fnfc;
      nfc = xf; fnfc = fx;
      xf = x; fx = fu,
      (* else *)
      If[x < xf, a = x, b = x];
      If[fu <= fnfc || nfc == xf,
        fulc = nfc; ffulc = fnfc;
        nfc = x; fnfc = fu,
        If[fu <= ffulc || fulc == xf || fulc == nfc,
          fulc = x; ffulc = fu;
        ]
      ]
    ];
    
    xm = 0.5*(a + b);
    tol1 = sqrtEps*Abs[xf] + xatol/3.0;
    tol2 = 2.0*tol1;
  ];
  
  {fx, xf}
];

(* ===================== Core Functions ===================== *)

(* SigmaSq: <h|h> = 4 deltaF Sum[|h[k]|^2 / Sn[k], {k, kmin, kmax-1}] *)
GWSigmaSq[hTilde_List, psd_List, deltaF_, fLow_, fHigh_] := Module[
  {nFreq = Length[hTilde], nTime, kmin, kmax, indices, mag2, psdSlice, result},
  nTime = (nFreq - 1) * 2;
  {kmin, kmax} = GetCutoffIndices[fLow, fHigh, deltaF, nTime];
  (* kmin..kmax-1 in 1-based indexing, matching pyCBC's kmin:kmax slice *)
  indices = Range[kmin, kmax - 1];
  mag2 = Abs[hTilde[[indices]]]^2;
  psdSlice = psd[[indices]];
  result = 4.0 * deltaF * Total[mag2 / psdSlice];
  Re[result] (* should be real, but ensure *)
];

GWSigma[hTilde_List, psd_List, deltaF_, fLow_, fHigh_] := 
  Sqrt[GWSigmaSq[hTilde, psd, deltaF, fLow, fHigh]];

(* Complex inner product: 4 deltaF Sum[Conj[h1[k]] * h2[k] / Sn[k]] *)
GWInnerProduct[h1Tilde_List, h2Tilde_List, psd_List, deltaF_, fLow_, fHigh_] := Module[
  {nFreq, nTime, kmin, kmax, indices, h1Slice, h2Slice, psdSlice},
  nFreq = Length[h1Tilde];
  nTime = (nFreq - 1) * 2;
  {kmin, kmax} = GetCutoffIndices[fLow, fHigh, deltaF, nTime];
  indices = Range[kmin, kmax - 1];
  h1Slice = h1Tilde[[indices]];
  h2Slice = h2Tilde[[indices]];
  psdSlice = psd[[indices]];
  4.0 * deltaF * Total[Conjugate[h1Slice] * h2Slice / psdSlice]
];

(* Normalized overlap (no time/phase maximization) *)
GWOverlap[h1Tilde_List, h2Tilde_List, psd_List, deltaF_, fLow_, fHigh_] := Module[
  {inner, sig1, sig2},
  inner = GWInnerProduct[h1Tilde, h2Tilde, psd, deltaF, fLow, fHigh];
  sig1 = GWSigma[h1Tilde, psd, deltaF, fLow, fHigh];
  sig2 = GWSigma[h2Tilde, psd, deltaF, fLow, fHigh];
  Re[inner] / (sig1 * sig2)
];

(* Cyclic time shift in frequency domain: h(f) -> h(f) * Exp[2 Pi I f dt] *)
CyclicTimeShift[hTilde_List, dt_, deltaF_] := Module[
  {nFreq = Length[hTilde], freqs, shift},
  freqs = Table[(k - 1) * deltaF, {k, 1, nFreq}];
  shift = Exp[2.0 Pi I freqs * dt];
  hTilde * shift
];

(* ===================== Match (Discrete) ===================== *)

(* Quadratic interpolation of peak, matching pyCBC's quadratic_interpolate_peak *)
QuadraticInterpolatePeak[left_, middle_, right_] := Module[
  {binOffset, peakValue},
  binOffset = 0.5 * (left - right) / (left - 2.0 * middle + right);
  peakValue = middle - 0.25 * (left - right) * binOffset;
  {binOffset, peakValue}
];

(* matched_filter_core equivalent:
   Computes IFFT of Conj[h1]*h2/PSD, returns the complex SNR time series
   and normalization factor. *)
MatchedFilterCore[h1Tilde_List, h2Tilde_List, psd_List, deltaF_, fLow_, fHigh_, 
                   h1Norm_: None] := Module[
  {nFreq, nTime, kmin, kmax, qtildeFull, snrTimeSeries, hNorm, norm, indices},
  nFreq = Length[h1Tilde];
  nTime = (nFreq - 1) * 2;
  {kmin, kmax} = GetCutoffIndices[fLow, fHigh, deltaF, nTime];
  
  (* Build the correlation in a full length-N complex array.
     Only positive frequency bins [kmin, kmax-1] are filled;
     negative frequency bins and out-of-band bins remain zero.
     This matches pyCBC's matched_filter_core exactly. *)
  qtildeFull = ConstantArray[0.0 + 0.0 I, nTime];
  indices = Range[kmin, kmax - 1];
  qtildeFull[[indices]] = Conjugate[h1Tilde[[indices]]] * h2Tilde[[indices]] / psd[[indices]];
  
  (* c2c IFFT: unnormalized inverse DFT, matching pyCBC's FFTW backward *)
  snrTimeSeries = InverseFFTc2c[qtildeFull];
  
  (* Normalization *)
  If[h1Norm === None,
    hNorm = GWSigmaSq[h1Tilde, psd, deltaF, fLow, fHigh],
    hNorm = h1Norm
  ];
  norm = (4.0 * deltaF) / Sqrt[hNorm];
  
  {snrTimeSeries, norm, hNorm}
];

(* Main match function: maximizes over time (via IFFT) and phase (via Abs) *)
Options[GWMatch] = {"SubsampleInterpolation" -> True};
GWMatch[h1In_List, h2In_List, deltaT_, psd_List, deltaF_, fLow_, fHigh_, 
        OptionsPattern[]] := Module[
  {h1 = h1In, h2 = h2In, nTime, h1Tilde, h2Tilde, snr, norm, hNorm, v2Norm, 
   maxSNR, maxID, matchVal, left, middle, right, idShift, subsample},
  
  If[Length[h1] != Length[h2],
    Message[GWMatch::lengthmismatch]; Return[$Failed]
  ];
  (* Zero-pad to even length if necessary *)
  If[OddQ[Length[h1]],
    h1 = Append[h1, 0.0];
    h2 = Append[h2, 0.0];
  ];
  nTime = Length[h1];
  
  (* FFT both waveforms *)
  h1Tilde = ForwardFFT[h1, deltaT];
  h2Tilde = ForwardFFT[h2, deltaT];
  
  (* Compute matched filter *)
  {snr, norm, hNorm} = MatchedFilterCore[h1Tilde, h2Tilde, psd, deltaF, fLow, fHigh];
  
  (* Find maximum of |SNR| *)
  (* maxID is the 1-based index of the peak *)
  maxID = First[Ordering[Abs[snr], -1]];
  maxSNR = Abs[snr[[maxID]]];
  
  (* v2 normalization *)
  v2Norm = GWSigmaSq[h2Tilde, psd, deltaF, fLow, fHigh];
  
  subsample = OptionValue["SubsampleInterpolation"];
  
  If[subsample,
    (* Quadratic interpolation of peak *)
    left = If[maxID == 1, Abs[snr[[-1]]], Abs[snr[[maxID - 1]]]];
    middle = maxSNR;
    right = If[maxID == nTime, Abs[snr[[1]]], Abs[snr[[maxID + 1]]]];
    {idShift, maxSNR} = QuadraticInterpolatePeak[left, middle, right];
    maxID = maxID + idShift;
  ];
  
  matchVal = maxSNR * norm / Sqrt[v2Norm];
  {matchVal, maxID}
];

(* ===================== Optimized Match ===================== *)

(* Optimized match: uses Brent's method for sub-sample accuracy.
   This is a faithful port of pyCBC's optimized_match function. *)
GWOptimizedMatch[h1In_List, h2In_List, deltaT_, psd_List, deltaF_, fLow_, fHigh_] := Module[
  {h1 = h1In, h2 = h2In, nTime, h1Tilde, h2Tilde, nFreq,
   (* First pass: discrete match to find approximate peak *)
   discreteMatch, discreteMaxID,
   (* After coarse time shift *)
   h2Shifted, h2ShiftedTilde,
   (* Frequency-domain quantities for optimization *)
   waveform1, waveform2, frequencies, psdArr,
   kmin, kmax, mask,
   (* Norms *)
   norm1, norm2, normProd,
   (* Optimization *)
   productFn, toMinimize, result, optDt, optMatch, optAngle},
  
  If[Length[h1] != Length[h2],
    Message[GWOptimizedMatch::lengthmismatch]; Return[$Failed]
  ];
  (* Zero-pad to even length if necessary *)
  If[OddQ[Length[h1]],
    h1 = Append[h1, 0.0];
    h2 = Append[h2, 0.0];
  ];
  nTime = Length[h1];
  
  (* FFT both waveforms *)
  h1Tilde = ForwardFFT[h1, deltaT];
  h2Tilde = ForwardFFT[h2, deltaT];
  nFreq = Length[h1Tilde];
  
  (* Step 1: Coarse match to find the best discrete time shift *)
  (* This replicates: _, max_id, _ = match(htilde, stilde, psd, ..., return_phase=True) *)
  Block[{snr, norm, hNorm, maxID, maxSNR, v2Norm},
    {snr, norm, hNorm} = MatchedFilterCore[h1Tilde, h2Tilde, psd, deltaF, fLow, fHigh];
    maxID = First[Ordering[Abs[snr], -1]];
    discreteMaxID = maxID;
  ];
  
  (* Step 2: Apply the coarse cyclic time shift to h2 *)
  (* pyCBC does: stilde = stilde.cyclic_time_shift(-max_id * delta_t)
     pyCBC's cyclic_time_shift(dt) multiplies by exp(-2j*pi*f*dt),
     so cyclic_time_shift(-max_id*dt) = exp(+2j*pi*f*max_id*dt).
     Our CyclicTimeShift multiplies by exp(+2j*pi*f*dt), so we pass
     +(maxID-1)*deltaT to get the same exp(+2j*pi*f*max_id_0*deltaT). *)
  h2ShiftedTilde = CyclicTimeShift[h2Tilde, +(discreteMaxID - 1) * deltaT, deltaF];
  
  (* Step 3: Extract frequency-domain data in the relevant band *)
  {kmin, kmax} = GetCutoffIndices[fLow, fHigh, deltaF, nTime];
  mask = Range[kmin, kmax - 1];
  
  waveform1 = h1Tilde[[mask]];
  waveform2 = h2ShiftedTilde[[mask]];
  frequencies = (mask - 1) * deltaF; (* 0-based frequencies *)
  psdArr = psd[[mask]];
  
  (* Step 4: Define the inner product as a function of time offset dt *)
  productFn[dt_?NumericQ] := Module[{offset, integral, mag, angle},
    offset = Exp[2.0 Pi I frequencies * dt];
    integral = Total[Conjugate[waveform1] * waveform2 * offset / psdArr] * deltaF;
    mag = 4.0 * Abs[integral];
    angle = Arg[integral];
    {mag, angle}
  ];
  
  toMinimize[dt_?NumericQ] := -First[productFn[dt]];
  
  (* Step 5: Compute normalizations *)
  norm1 = GWSigmaSq[h1Tilde, psd, deltaF, fLow, fHigh];
  norm2 = GWSigmaSq[h2Tilde, psd, deltaF, fLow, fHigh]; (* shift doesn't change norm *)
  normProd = Sqrt[norm1 * norm2];
  
  (* Step 6: Optimize over dt in (-deltaT, deltaT) using Brent's method *)
  (* This exactly matches pyCBC's use of scipy.optimize.minimize_scalar(method='bounded') *)
  result = BrentMinimize[toMinimize, {-deltaT, deltaT}];
  
  optDt = result[[2]];
  optMatch = -result[[1]] / normProd;
  
  {optMatch, optDt / deltaT + (discreteMaxID - 1)}
];

(* ===================== Mismatch ===================== *)

GWMismatch[h1_List, h2_List, deltaT_, psd_List, deltaF_, fLow_, fHigh_] := Module[
  {result},
  result = GWOptimizedMatch[h1, h2, deltaT, psd, deltaF, fLow, fHigh];
  {1.0 - result[[1]], result[[2]]}
];

(* ===================== PSD Loading ===================== *)

(* Load a two-column (freq, PSD) text file and interpolate onto a regular grid *)
LoadPSD[filename_String, deltaF_, fLow_, nFreq_Integer] := Module[
  {data, interpFn, freqs, psdValues, kmin},
  
  (* Read two-column data *)
  data = Import[filename, "Table"];
  If[data === $Failed, Message[LoadPSD::filenotfound]; Return[$Failed]];
  
  (* Create interpolation function (log-log for PSD, as in LALSim) *)
  interpFn = Interpolation[
    Transpose[{Log[data[[All, 1]]], Log[data[[All, 2]]]}],
    InterpolationOrder -> 1
  ];
  
  (* Evaluate on the target frequency grid *)
  freqs = Table[(k - 1) * deltaF, {k, 1, nFreq}];
  psdValues = Table[
    If[f <= 0 || f < fLow,
      Infinity,
      If[f >= data[[-1, 1]], (* beyond data range *)
        Infinity,
        Exp[interpFn[Log[f]]]
      ]
    ],
    {f, freqs}
  ];
  
  psdValues
];

(* Flat (unity) PSD *)
FlatPSD[nFreq_Integer, deltaF_, fLow_] := Module[
  {psd, kmin},
  psd = ConstantArray[1.0, nFreq];
  kmin = Max[1, Floor[fLow / deltaF]];
  psd[[1 ;; kmin]] = Infinity;
  psd
];

(* ===================== Detector Noise PSD Models ===================== *)

(* Physical constants *)
$cSI = 299792458.0;            (* speed of light [m/s] *)
$hbarSI = 1.054571817*^-34;    (* reduced Planck constant [J s] *)
$kBSI = 1.380649*^-23;         (* Boltzmann constant [J/K] *)

(* --- Component noise models (from LALSimNoisePSD.c) --- *)

(* Suspension thermal noise: broadband above f0.
   S_h(f) = 2 k T / (L^2 M Q (Pi f0)^3) * (f0/f)^5 *)
SuspThermalPSD[f_, L_, M_, T_, f0_, Q_] :=
  2.0 $kBSI T / (L^2 M Q (Pi f0)^3) (f0 / f)^5;

(* Coating/mirror thermal noise: broadband below resonance.
   S_h(f) = 2 k T / (L^2 M Q (Pi f0)^3) * (f0/f) *)
CoatThermalPSD[f_, L_, M_, T_, f0_, Q_] :=
  2.0 $kBSI T / (L^2 M Q (Pi f0)^3) (f0 / f);

(* aLIGO thermal noise: suspension + coating *)
aLIGOThermalPSD[f_] := Module[{L, M, T, fS, QS, fC, QC},
  L = 3995.0; M = 40.0; T = 290.0;
  fS = 9.0; QS = 6.0*^10; fC = 1.0*^4; QC = 6.0*^6;
  SuspThermalPSD[f, L, M, T, fS, QS] + CoatThermalPSD[f, L, M, T, fC, QC]
];

(* Quantum noise: Buonanno & Chen, Phys. Rev. D 64, 042006 (2001).
   Full radiation-pressure + shot noise model.
   Adapted from GWINC / LALSimNoisePSD.c XLALSimNoisePSDQuantum. *)
QuantumNoisePSD[f_?NumericQ, i0_, lambda_, L_, M_, mirA_, bsA_, tITM_, tPRM_, tSRM_,
                ds_, zeta_, eta_] :=
 Module[
  {omega, omega0, lSR, lPD, tau, rho, phi, lArm, gAC, eps, r1, rarm,
   gPRC, iCirc, iSQL, kp, bt, hSQL,
   c11, c12, c21, d1, d2, p11, p12, p21, q11, n11, n12, n21, n22},

  omega  = 2.0 Pi f;
  omega0 = 2.0 Pi $cSI / lambda;
  lSR    = bsA;
  lPD    = 1.0 - eta;
  tau    = Sqrt[tSRM];
  rho    = Sqrt[1.0 - tSRM];
  phi    = (Pi - ds) / 2.0;
  lArm   = 2.0 mirA;
  gAC    = tITM $cSI / (4.0 L);
  eps    = lArm / (2.0 gAC L / $cSI);
  r1     = Sqrt[1.0 - tITM];
  rarm   = r1 - tITM Sqrt[1.0 - 2.0 mirA] / (1.0 - r1 Sqrt[1.0 - 2.0 mirA]);
  gPRC   = tPRM / (1.0 + Sqrt[1.0 - tPRM] rarm Sqrt[1.0 - bsA])^2;
  iCirc  = gPRC i0;
  iSQL   = M L^2 gAC^4 / (4.0 omega0);
  kp     = 2.0 (iCirc / iSQL) gAC^4 / (omega^2 (gAC^2 + omega^2));
  bt     = ArcTan[omega / gAC];
  hSQL   = Sqrt[8.0 $hbarSI / (M (omega L)^2)];

  (* Transfer-matrix elements [BC Eqs. 5.8-5.12] *)
  c11 = Sqrt[1 - lPD] (
    (1 + rho^2) (Cos[2 phi] + kp/2 Sin[2 phi]) - 2 rho Cos[2 bt] -
    eps/4 (-2 (1 + Exp[2 I bt])^2 rho +
      4 (1 + rho^2) Cos[bt]^2 Cos[2 phi] +
      (3 + Exp[2 I bt]) kp (1 + rho^2) Sin[2 phi]) +
    lSR (Exp[2 I bt] rho -
      (1 + rho^2)/2 (Cos[2 phi] + kp/2 Sin[2 phi]))
  );
  (* c22 = c11 *)

  c12 = Sqrt[1 - lPD] tau^2 (
    -(Sin[2 phi] + kp Sin[phi]^2) +
    eps/2 Sin[phi] ((3 + Exp[2 I bt]) kp Sin[phi] +
      4 Cos[bt]^2 Cos[phi]) +
    lSR/2 (Sin[2 phi] + kp Sin[phi]^2)
  );

  c21 = Sqrt[1 - lPD] tau^2 (
    (Sin[2 phi] - kp Cos[phi]^2) +
    eps/2 Cos[phi] ((3 + Exp[2 I bt]) kp Sin[phi] -
      4 Cos[bt]^2 Sin[phi]) +
    lSR/2 (-Sin[2 phi] + kp Cos[phi]^2)
  );

  d1 = Sqrt[1 - lPD] (
    -(1 + rho Exp[2 I bt]) Sin[phi] +
    eps/4 (3 + rho + 2 rho Exp[4 I bt] +
      Exp[2 I bt] (1 + 5 rho)) Sin[phi] +
    lSR/2 Exp[2 I bt] rho Sin[phi]
  );

  d2 = Sqrt[1 - lPD] (
    -(-1 + rho Exp[2 I bt]) Cos[phi] +
    eps/4 (-3 + rho + 2 rho Exp[4 I bt] +
      Exp[2 I bt] (-1 + 5 rho)) Cos[phi] +
    lSR/2 Exp[2 I bt] rho Cos[phi]
  );

  p11 = Sqrt[1 - lPD] Sqrt[lSR] tau / 2 (
    -2 rho Exp[2 I bt] + 2 Cos[2 phi] + kp Sin[2 phi]
  );
  (* p22 = p11 *)

  p12 = -Sqrt[1 - lPD] Sqrt[lSR] tau Sin[phi] (
    2 Cos[phi] + kp Sin[phi]
  );

  p21 = Sqrt[1 - lPD] Sqrt[lSR] tau Cos[phi] (
    2 Sin[phi] - kp Cos[phi]
  );

  q11 = Sqrt[lPD] (
    Exp[-2 I bt] + rho^2 Exp[2 I bt] -
    rho (2 Cos[2 phi] + kp Sin[2 phi]) +
    eps rho / 2 (
      Exp[-2 I bt] Cos[2 phi] +
      Exp[2 I bt] (-2 rho - 2 rho Cos[2 bt] +
        Cos[2 phi] + kp Sin[2 phi]) +
      2 Cos[2 phi] + 3 kp Sin[2 phi]) -
    lSR rho / 2 (2 rho Exp[2 I bt] -
      2 Cos[2 phi] - kp Sin[2 phi])
  );
  (* q22 = q11, q12 = q21 = 0 *)

  n11 = Sqrt[1 - lPD] Sqrt[eps / 2] tau (
    kp (1 + rho Exp[2 I bt]) Sin[phi] +
    2 Cos[bt] (Exp[-I bt] Cos[phi] -
      rho Exp[I bt] (Cos[phi] + kp Sin[phi]))
  );

  n22 = -Sqrt[1 - lPD] Sqrt[2 eps] tau (
    -Exp[-I bt] + rho Exp[I bt]
  ) Cos[bt] Cos[phi];

  n12 = -Sqrt[1 - lPD] Sqrt[2 eps] tau (
    Exp[-I bt] + rho Exp[I bt]
  ) Cos[bt] Sin[phi];

  n21 = Sqrt[1 - lPD] Sqrt[eps / 2] tau (
    -kp (1 + rho) Cos[phi] +
    2 Cos[bt] (Exp[-I bt] + rho Exp[I bt]) Cos[bt] Sin[phi]
  );

  (* PSD [BC Eq. 5.13] *)
  Re[hSQL^2 / (2 kp tau^2 Abs[d1 Sin[zeta] + d2 Cos[zeta]]^2) (
    Abs[c11 Sin[zeta] + c21 Cos[zeta]]^2 +
    Abs[c12 Sin[zeta] + c11 Cos[zeta]]^2 +
    Abs[p11 Sin[zeta] + p21 Cos[zeta]]^2 +
    Abs[p12 Sin[zeta] + p11 Cos[zeta]]^2 +
    Abs[q11]^2 +
    Abs[n11 Sin[zeta] + n21 Cos[zeta]]^2 +
    Abs[n12 Sin[zeta] + n22 Cos[zeta]]^2
  )]
];

(* --- aLIGO configurations (from LIGO-T0900288-v3 / LIGO-T070247-01) ---
   All use thermal + quantum noise; valid above ~9 Hz.
   Parameters: L=3995m, M=40kg, lambda=1.064um, A=37.5e-6, A_BS=0.002,
   T_ITM=0.014, T_PRM=0.027, T_SRM=0.2, eta=0.9 *)

aLIGOQuantumPSD[f_?NumericQ, i0_, ds_, zeta_] :=
  QuantumNoisePSD[f, i0, 1.064*^-6, 3995.0, 40.0, 37.5*^-6, 0.002,
    0.014, 0.027, 0.2, ds, zeta, 0.9];

(* Zero Detuned, High Power (125 W) - THE standard aLIGO design sensitivity *)
PSDaLIGOZeroDetHighPower[f_?NumericQ] :=
  aLIGOThermalPSD[f] + aLIGOQuantumPSD[f, 125.0, 0.0, 116.0 Degree];

(* Zero Detuned, Low Power (25 W) *)
PSDaLIGOZeroDetLowPower[f_?NumericQ] :=
  aLIGOThermalPSD[f] + aLIGOQuantumPSD[f, 25.0, 0.0, 116.0 Degree];

(* NS-NS Optimized *)
PSDaLIGONSNSOpt[f_?NumericQ] :=
  aLIGOThermalPSD[f] + aLIGOQuantumPSD[f, 125.0, 11.0 Degree, 103.0 Degree];

(* BHBH 20-degree Detune (20 W input) *)
PSDaLIGOBHBH20Deg[f_?NumericQ] :=
  aLIGOThermalPSD[f] +
  QuantumNoisePSD[f, 20.0, 1.064*^-6, 3995.0, 40.0, 37.5*^-6, 0.002,
    0.014, 0.027, 0.2, 20.0 Degree, 105.0 Degree, 0.9];

(* High Frequency (T_SRM=0.011) *)
PSDaLIGOHighFrequency[f_?NumericQ] :=
  aLIGOThermalPSD[f] +
  QuantumNoisePSD[f, 125.0, 1.064*^-6, 3995.0, 40.0, 37.5*^-6, 0.002,
    0.014, 0.027, 0.011, 4.7 Degree, 128.0 Degree, 0.9];

(* --- Phenomenological PSD fits --- *)

(* Advanced Virgo design sensitivity.
   Phenomenological fit from Miao, Zhao, Chen, PRD 86, 104003 (2012).
   Reference: http://wwwcascina.virgo.infn.it/advirgo *)
PSDAdvVirgo[f_?NumericQ] := Module[{x, x2, asd},
  x = Log[f / 300.0];
  x2 = x^2;
  asd = 1.259*^-24 (
    0.07 Exp[-0.142 - 1.437 x + 0.407 x2] +
    3.1  Exp[-0.466 - 1.043 x - 0.548 x2] +
    0.4  Exp[-0.304 + 2.896 x - 0.293 x2] +
    0.09 Exp[ 1.466 + 3.722 x - 0.984 x2]);
  asd^2
];

(* KAGRA design sensitivity.
   Phenomenological fit from Miao, Zhao, Chen, PRD 86, 104003 (2012).
   Reference: http://gwcenter.icrr.u-tokyo.ac.jp/en/researcher/parameter *)
PSDKAGRA[f_?NumericQ] := Module[{x, x2, asd},
  x = Log[f / 100.0];
  x2 = x^2;
  asd = 6.499*^-25 (
    9.72*^-9 Exp[-1.43 - 9.88 x - 0.23 x2] +
    1.17     Exp[ 0.14 - 3.10 x - 0.26 x2] +
    1.70     Exp[ 0.14 + 1.09 x - 0.013 x2] +
    1.25     Exp[ 0.071 + 2.83 x - 4.91 x2]);
  asd^2
];

(* Initial Virgo design sensitivity.
   Phenomenological fit from the Virgo website.
   Reference: LALSimNoisePSD.c XLALSimNoisePSDVirgo *)
PSDVirgo[f_?NumericQ] := Module[{x, s0},
  x = f / 500.0;
  s0 = 10.2*^-46;
  s0 ((7.87 x)^(-4.8) + 6.0/17.0/x + 1.0 + x^2)
];

(* Initial LIGO Science Requirements Document (SRD) sensitivity.
   Phenomenological fit to the SRD curve.
   Reference: LALSimNoisePSD.c XLALSimNoisePSDiLIGOSRD *)
PSDiLIGOSRD[f_?NumericQ] := Module[
  {aseis = 1.57271, pseis = -14.0,
   athrm = 3.80591*^-19, pthrm = -2.0,
   ashot = 1.12277*^-23, fshot = 89.3676},
  aseis^2 f^(2 pseis) + athrm^2 f^(2 pthrm) +
    ashot^2 (1.0 + (f / fshot)^2)
];

(* ================================================================== *)
(* LISA PSD model                                                      *)
(* Based on Robson, Cornish, Liu, Class. Quant. Grav. 36 (2019) 105011*)
(* arXiv:1803.01944                                                    *)
(* Code: github.com/eXtremeGravityInstitute/LISA_Sensitivity           *)
(* ================================================================== *)

(* LISA design parameters *)
$LISAArmLength = 2.5*^9;                                   (* arm length [m] *)
$LISANumChannels = 2;                                       (* number of TDI channels *)
$LISAfstar = $cSI / (2.0 Pi $LISAArmLength);               (* transfer frequency ~19.09 mHz *)

(* Single-link optical metrology noise, Eq. (10) *)
lisaPoms[f_] := (1.5*^-11)^2 (1.0 + (2.0*^-3 / f)^4);

(* Single test mass acceleration noise, Eq. (11) *)
lisaPacc[f_] := (3.0*^-15)^2 (1.0 + (0.4*^-3 / f)^2) (1.0 + (f / 8.0*^-3)^4);

(* Strain PSD in a Michelson-style TDI channel, Eq. (12) *)
lisaPn[f_] :=
  (lisaPoms[f] + 2.0 (1.0 + Cos[f / $LISAfstar]^2) lisaPacc[f] / (2.0 Pi f)^4) /
    $LISAArmLength^2;

(* Approximate sky- and polarization-averaged response function, Eq. (9) *)
lisaR[f_] :=
  (3.0 / 20.0) / (1.0 + 0.6 (f / $LISAfstar)^2) $LISANumChannels;

(* Galactic binary confusion noise estimate, Eq. (14) and Table I.
   tObsYears is the observation time in years. *)
lisaSnC[f_, tObsYears_] := Module[
  {alpha, beta, kappa, gamma, fKnee, a, sc},
  {alpha, beta, kappa, gamma, fKnee} = Which[
    tObsYears < 0.75, {0.133,  243., 482.,  917., 2.58*^-3},
    tObsYears < 1.5,  {0.171,  292., 1020., 1680., 2.15*^-3},
    tObsYears < 3.0,  {0.165,  299., 611.,  1340., 1.73*^-3},
    True,             {0.138, -221., 521.,  1680., 1.13*^-3}
  ];
  a = 1.8*^-44 / $LISANumChannels;
  sc = (1.0 + Tanh[gamma (fKnee - f)]) *
       Exp[-f^alpha + beta f Sin[kappa f]] a f^(-7.0/3.0);
  sc
];

(* LISA sky-averaged sensitivity (no confusion noise): Sn = Pn/R *)
PSDLISA[f_?NumericQ] := lisaPn[f] / lisaR[f];

(* LISA sensitivity with galactic confusion noise for different mission durations *)
PSDLISAConfusion05yr[f_?NumericQ] := lisaPn[f] / lisaR[f] + lisaSnC[f, 0.5];
PSDLISAConfusion1yr[f_?NumericQ]  := lisaPn[f] / lisaR[f] + lisaSnC[f, 1.0];
PSDLISAConfusion2yr[f_?NumericQ]  := lisaPn[f] / lisaR[f] + lisaSnC[f, 2.0];
PSDLISAConfusion4yr[f_?NumericQ]  := lisaPn[f] / lisaR[f] + lisaSnC[f, 4.0];

(* --- Tabulated PSD models (from LIGO-P1600143-v18 data files) --- *)

(* Directory containing tabulated PSD data files.
   Defaults to psd_data/ subdirectory relative to this package file. *)
$GWMMatchPSDDataDirectory = FileNameJoin[{DirectoryName[$InputFileName], "psd_data"}];

(* Cache for loaded interpolation functions *)
$tabulatedPSDCache = <||>;

(* Load a two-column (frequency, ASD) file and return a PSD function.
   Uses log-log cubic spline interpolation, matching LALSimNoisePSD.c.
   The ASD values are squared to give PSD. Results are cached. *)
loadTabulatedASD[filename_String] := Module[
  {filepath, data, logFreqs, logASDs, interpFn, fMin, fMax, psdFn},

  (* Return cached version if available *)
  If[KeyExistsQ[$tabulatedPSDCache, filename],
    Return[$tabulatedPSDCache[filename]]
  ];

  filepath = FileNameJoin[{$GWMMatchPSDDataDirectory, filename}];
  data = Import[filepath, "Table"];
  If[data === $Failed || data === {},
    Message[LoadPSD::filenotfound, filepath];
    Return[$Failed]
  ];

  (* Build log-log cubic spline interpolation (matching LALSim gsl_interp_cspline) *)
  logFreqs = Log[data[[All, 1]]];
  logASDs = Log[data[[All, 2]]];
  interpFn = Interpolation[
    Transpose[{logFreqs, logASDs}],
    InterpolationOrder -> 3
  ];

  fMin = data[[1, 1]];
  fMax = data[[-1, 1]];

  (* Return a function: f -> PSD(f) = ASD(f)^2 *)
  psdFn = Function[f,
    If[f < fMin || f > fMax,
      Infinity,
      Exp[2.0 interpFn[Log[f]]]  (* ASD^2 = Exp[2 * log(ASD)] *)
    ]
  ];

  (* Cache and return *)
  $tabulatedPSDCache[filename] = psdFn;
  psdFn
];

(* Einstein Telescope D configuration sensitivity.
   Reference: LIGO-P1600143-v18, Hild et al. Class. Quantum Grav. 28 (2011) 094013 *)
PSDEinsteinTelescopeP1600143[f_?NumericQ] :=
  loadTabulatedASD["LIGO-P1600143-v18-ET_D.txt"][f];

(* Cosmic Explorer baseline sensitivity.
   Reference: LIGO-P1600143-v18 *)
PSDCosmicExplorerP1600143[f_?NumericQ] :=
  loadTabulatedASD["LIGO-P1600143-v18-CE.txt"][f];

(* Cosmic Explorer pessimistic sensitivity.
   Reference: LIGO-P1600143-v18 *)
PSDCosmicExplorerPessimisticP1600143[f_?NumericQ] :=
  loadTabulatedASD["LIGO-P1600143-v18-CE_Pessimistic.txt"][f];

(* Cosmic Explorer wideband sensitivity.
   Reference: LIGO-P1600143-v18 *)
PSDCosmicExplorerWidebandP1600143[f_?NumericQ] :=
  loadTabulatedASD["LIGO-P1600143-v18-CE_Wideband.txt"][f];

(* aLIGO zero-detuning high power GWINC (tabulated).
   This is the PSD used by pyCBC's aLIGOZeroDetHighPowerGWINC.
   Reference: LIGO-T0900288-v3 *)
PSDaLIGOZeroDetHighPowerGWINC[f_?NumericQ] :=
  loadTabulatedASD["LIGO-T0900288-v3-ZERO_DET_high_P.txt"][f];

(* --- Main interface --- *)

$detectorPSDFunctions = <|
  "aLIGOZeroDetHighPower" -> PSDaLIGOZeroDetHighPower,
  "aLIGOZeroDetLowPower"  -> PSDaLIGOZeroDetLowPower,
  "aLIGONSNSOpt"          -> PSDaLIGONSNSOpt,
  "aLIGOBHBH20Deg"        -> PSDaLIGOBHBH20Deg,
  "aLIGOHighFrequency"    -> PSDaLIGOHighFrequency,
  "AdvVirgo"              -> PSDAdvVirgo,
  "KAGRA"                 -> PSDKAGRA,
  "Virgo"                 -> PSDVirgo,
  "iLIGOSRD"              -> PSDiLIGOSRD,
  "EinsteinTelescopeP1600143"         -> PSDEinsteinTelescopeP1600143,
  "CosmicExplorerP1600143"            -> PSDCosmicExplorerP1600143,
  "CosmicExplorerPessimisticP1600143" -> PSDCosmicExplorerPessimisticP1600143,
  "CosmicExplorerWidebandP1600143"    -> PSDCosmicExplorerWidebandP1600143,
  "aLIGOZeroDetHighPowerGWINC"         -> PSDaLIGOZeroDetHighPowerGWINC,
  "LISA"                                -> PSDLISA,
  "LISAConfusion05yr"                   -> PSDLISAConfusion05yr,
  "LISAConfusion1yr"                    -> PSDLISAConfusion1yr,
  "LISAConfusion2yr"                    -> PSDLISAConfusion2yr,
  "LISAConfusion4yr"                    -> PSDLISAConfusion4yr
|>;

ListDetectorPSDs[] := Keys[$detectorPSDFunctions];

DetectorPSD::unknown = "Unknown detector PSD model `1`. Use ListDetectorPSDs[] for available models.";

DetectorPSD[name_String, nFreq_Integer, deltaF_, fLow_] := Module[
  {psdFunc, freqs, psd, kmin},

  If[!KeyExistsQ[$detectorPSDFunctions, name],
    Message[DetectorPSD::unknown, name]; Return[$Failed]
  ];
  psdFunc = $detectorPSDFunctions[name];

  freqs = Table[(k - 1) deltaF, {k, 1, nFreq}];
  kmin = Max[1, Ceiling[fLow / deltaF]];

  psd = Table[
    If[k <= kmin || freqs[[k]] <= 0,
      Infinity,
      psdFunc[freqs[[k]]]
    ],
    {k, 1, nFreq}
  ];
  psd
];

End[];
EndPackage[];
