(* =========================================================== *)
(*  01-BasicUsage.wls                                          *)
(*  GWMMatch — end-to-end match computation                   *)
(*                                                             *)
(*  Run:  wolframscript -file 01-BasicUsage.wls               *)
(*  Or:   Get["path/to/notebooks/01-BasicUsage.wls"]          *)
(* =========================================================== *)

<< GWMMatch`

(* ----------------------------------------------------------- *)
(* 1. Build two simple chirp-like waveforms                    *)
(* ----------------------------------------------------------- *)

sampleRate = 2048.0;   (* Hz *)
duration   = 2.0;      (* seconds *)
f0         = 40.0;     (* starting frequency [Hz] *)
fDot       = 30.0;     (* frequency sweep rate [Hz/s] *)

nTime  = Round[sampleRate * duration];
deltaT = 1.0 / sampleRate;
deltaF = 1.0 / duration;
fLow   = 20.0;   (* lower frequency cutoff [Hz] *)
fHigh  = 512.0;  (* upper frequency cutoff [Hz] *)

(* Reference waveform: a frequency-swept sinusoid with a Planck taper *)
tVals = Table[n * deltaT, {n, 0, nTime - 1}];
h1Raw = Table[Sin[2 Pi (f0 + fDot t) t], {t, tVals}];

(* Template: same waveform but shifted in time by 20 ms and scaled by 2 *)
dtShift      = 0.020;  (* seconds *)
shiftSamples = Round[dtShift * sampleRate];
h2Raw        = 2.0 * RotateLeft[h1Raw, shiftSamples];

(* Apply a Planck taper to suppress spectral leakage at the edges *)
h1 = GWTaper[h1Raw];
h2 = GWTaper[h2Raw];

Print["Waveform length: ", nTime, " samples  (", duration, " s at ", sampleRate, " Hz)"];
Print["Time shift applied to h2: ", shiftSamples, " samples (", dtShift * 1000, " ms)"];

(* ----------------------------------------------------------- *)
(* 2. Build a flat PSD (whitened data)                         *)
(* ----------------------------------------------------------- *)

nFreq  = Floor[nTime / 2] + 1;
psdFlat = FlatPSD[nFreq, deltaF, fLow];

(* ----------------------------------------------------------- *)
(* 3. Basic overlap (no time/phase maximization)               *)
(* ----------------------------------------------------------- *)

(* FFT the tapered waveforms *)
h1Tilde = GWMMatch`Private`ForwardFFT[h1, deltaT];
h2Tilde = GWMMatch`Private`ForwardFFT[h2, deltaT];

overlap = GWOverlap[h1Tilde, h2Tilde, psdFlat, deltaF, fLow, fHigh];
Print["\n-- GWOverlap (no maximization) --"];
Print["Overlap = ", N[overlap, 6]];
Print["(Expected: << 1 because of the time shift)"];

(* ----------------------------------------------------------- *)
(* 4. Discrete match (maximized over time and phase)           *)
(* ----------------------------------------------------------- *)

{matchDiscrete, peakIndex} = GWMatch[h1, h2, deltaT, psdFlat, deltaF, fLow, fHigh];
recoveredShift             = (peakIndex - 1) * deltaT;  (* seconds *)

Print["\n-- GWMatch (discrete) --"];
Print["Match       = ", N[matchDiscrete, 6]];
Print["Peak index  = ", N[peakIndex, 8]];
Print["Recovered time shift = ", N[recoveredShift * 1000, 6], " ms  (true: ", dtShift * 1000, " ms)"];

(* ----------------------------------------------------------- *)
(* 5. Optimized match (sub-sample accuracy via Brent's method) *)
(* ----------------------------------------------------------- *)

{matchOpt, peakIndexOpt} = GWOptimizedMatch[h1, h2, deltaT, psdFlat, deltaF, fLow, fHigh];
recoveredShiftOpt         = (peakIndexOpt - 1) * deltaT;

Print["\n-- GWOptimizedMatch (Brent sub-sample refinement) --"];
Print["Match       = ", N[matchOpt, 6]];
Print["Peak index  = ", N[peakIndexOpt, 8]];
Print["Recovered time shift = ", N[recoveredShiftOpt * 1000, 6], " ms  (true: ", dtShift * 1000, " ms)"];

(* ----------------------------------------------------------- *)
(* 6. Mismatch                                                 *)
(* ----------------------------------------------------------- *)

{mm, _} = GWMismatch[h1, h2, deltaT, psdFlat, deltaF, fLow, fHigh];
Print["\n-- GWMismatch --"];
Print["Mismatch = ", N[mm, 6], "  (= 1 - match)"];

(* ----------------------------------------------------------- *)
(* 7. Repeat with a realistic aLIGO PSD                       *)
(* ----------------------------------------------------------- *)

psdALIGO = DetectorPSD["aLIGOZeroDetHighPower", nFreq, deltaF, fLow];

{matchALIGO, _} = GWOptimizedMatch[h1, h2, deltaT, psdALIGO, deltaF, fLow, fHigh];
Print["\n-- GWOptimizedMatch with aLIGO PSD --"];
Print["Match = ", N[matchALIGO, 6]];
Print["Mismatch = ", N[1 - matchALIGO, 6]];

(* ----------------------------------------------------------- *)
(* 8. Self-match sanity check                                  *)
(* ----------------------------------------------------------- *)

{selfMatch, _} = GWOptimizedMatch[h1, h1, deltaT, psdALIGO, deltaF, fLow, fHigh];
Print["\n-- Self-match sanity check --"];
Print["match(h1, h1) = ", N[selfMatch, 10], "  (should be 1.0)"];
