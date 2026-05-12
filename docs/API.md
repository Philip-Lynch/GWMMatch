# GWMMatch API Reference

Full reference for all public functions in `GWMMatch`. Load the package with `<< GWMMatch``.

---

## Conventions

- All waveforms are **real-valued** lists sampled at interval `deltaT` seconds.
- Frequency-domain arrays (`hTilde`) are **one-sided** complex lists of length `Floor[nTime/2] + 1`, indexed from bin 0 (DC) to bin Nyquist.
- PSD arrays have the same length as frequency-domain arrays; bins below `fLow` and the DC bin are set to `Infinity` to exclude them from inner products.
- Indices are **1-based** throughout (Mathematica convention).
- The match is normalised to `[0, 1]`; the mismatch is `1 - match`.

---

## Inner Products and Norms

### `GWInnerProduct`

```mathematica
GWInnerProduct[h1Tilde, h2Tilde, psd, deltaF, fLow, fHigh]
```

Computes the complex noise-weighted inner product:

$$\langle h_1 | h_2 \rangle = 4\,\Delta f \sum_{k=k_\min}^{k_\max - 1} \frac{\tilde{h}_1^*(f_k)\,\tilde{h}_2(f_k)}{S_n(f_k)}$$

**Arguments**

| Name | Type | Description |
|---|---|---|
| `h1Tilde` | complex list | One-sided FFT of waveform 1 |
| `h2Tilde` | complex list | One-sided FFT of waveform 2 |
| `psd` | real list | One-sided PSD, same length as `h1Tilde`; `Infinity` outside band |
| `deltaF` | real | Frequency bin spacing [Hz] |
| `fLow` | real | Lower frequency cutoff [Hz] |
| `fHigh` | real | Upper frequency cutoff [Hz] |

**Returns** a complex number. Take `Re[...]` for the real inner product.

---

### `GWSigmaSq`

```mathematica
GWSigmaSq[hTilde, psd, deltaF, fLow, fHigh]
```

Signal norm squared: $\sigma^2 = \langle h | h \rangle$. Returns a non-negative real number.

---

### `GWSigma`

```mathematica
GWSigma[hTilde, psd, deltaF, fLow, fHigh]
```

Signal norm: $\sigma = \sqrt{\langle h | h \rangle}$.

---

### `GWOverlap`

```mathematica
GWOverlap[h1Tilde, h2Tilde, psd, deltaF, fLow, fHigh]
```

Normalized overlap without time or phase maximization:

$$\mathcal{O}(h_1, h_2) = \frac{\mathrm{Re}\langle h_1 | h_2 \rangle}{\sigma_1 \sigma_2} \in [-1, 1]$$

Takes **frequency-domain** inputs. Use `GWMatch` or `GWOptimizedMatch` for time-domain inputs with maximization.

---

## Match Functions

All match functions take **time-domain** waveforms as input.

### `GWMatch`

```mathematica
GWMatch[h1, h2, deltaT, psd, deltaF, fLow, fHigh]
GWMatch[h1, h2, deltaT, psd, deltaF, fLow, fHigh, "SubsampleInterpolation" -> True]
```

Computes the match maximized over **discrete** time shifts (via an IFFT of the cross-correlation) and all phase shifts (via `Abs`). Optionally applies quadratic sub-sample interpolation at the peak.

**Arguments**

| Name | Type | Description |
|---|---|---|
| `h1`, `h2` | real lists | Time-domain waveforms, must have equal length |
| `deltaT` | real | Sample spacing [s] |
| `psd` | real list | One-sided PSD of length `Floor[Length[h1]/2] + 1` |
| `deltaF` | real | Frequency bin spacing [Hz] (= `1 / (Length[h1] * deltaT)`) |
| `fLow` | real | Lower frequency cutoff [Hz] |
| `fHigh` | real | Upper frequency cutoff [Hz] |

**Options**

| Option | Default | Description |
|---|---|---|
| `"SubsampleInterpolation"` | `True` | Apply quadratic interpolation at the SNR peak |

**Returns** `{match, peakIndex}` where `match ∈ [0, 1]` and `peakIndex` is the 1-based sample index of the best time shift (fractional if sub-sample interpolation is on).

**Notes**
- Waveforms are zero-padded to even length if necessary.
- The phase is maximized analytically by taking `Abs` of the complex SNR time series.
- For higher accuracy use `GWOptimizedMatch`.

---

### `GWOptimizedMatch`

```mathematica
GWOptimizedMatch[h1, h2, deltaT, psd, deltaF, fLow, fHigh]
```

Computes the match using a two-step procedure:

1. **Coarse**: finds the best discrete time-shift bin with `GWMatch`.
2. **Fine**: applies Brent's bounded scalar minimizer within ±1 sample of the coarse peak to find the optimal continuous time shift.

This is a faithful port of pyCBC's `optimized_match`. It is more accurate than `GWMatch` at the cost of additional function evaluations.

**Returns** `{match, peakIndex}` where `peakIndex` is a real-valued sample index.

---

### `GWMismatch`

```mathematica
GWMismatch[h1, h2, deltaT, psd, deltaF, fLow, fHigh]
```

Convenience wrapper. Returns `{1 - match, peakIndex}` using `GWOptimizedMatch`.

---

## Tapering

### `GWTaper`

```mathematica
GWTaper[h]
GWTaper[h, method]
GWTaper[h, method, taperSamples]
```

Applies a **Planck-window taper** to suppress spectral leakage at the waveform boundaries. Ports LALSuite's `XLALSimInspiralREAL8WaveTaper` (McKechan, Robinson & Sathyaprakash, CQG 27 (2010) 084020).

**Arguments**

| Name | Type | Description |
|---|---|---|
| `h` | real list | Time-domain waveform |
| `method` | string | `"TAPER_START"`, `"TAPER_END"`, or `"TAPER_STARTEND"` (default) |
| `taperSamples` | integer, list, or `Automatic` | Taper width in samples |

**`taperSamples` forms**

| Value | Effect |
|---|---|
| `Automatic` (default) | Auto-detect width from second local extremum of `\|h\|`, matching LALSuite |
| `n` (integer) | Fixed width `n` at both ends |
| `{nStart, nEnd}` | Independent widths; either element can be `Automatic` |

**Returns** the tapered waveform (same length as input). The first and last nonzero samples are set to zero; the transition follows the Planck sigmoid.

**Example**

```mathematica
(* Auto-detect both ends *)
hTapered = GWTaper[h]

(* Fixed 200-sample taper at both ends *)
hTapered = GWTaper[h, "TAPER_STARTEND", 200]

(* Auto-detect start, fixed 150 samples at end *)
hTapered = GWTaper[h, "TAPER_STARTEND", {Automatic, 150}]

(* Taper the start only *)
hTapered = GWTaper[h, "TAPER_START"]
```

---

## PSD Utilities

### `DetectorPSD`

```mathematica
DetectorPSD[name, nFreq, deltaF, fLow]
```

Returns a one-sided PSD array of length `nFreq` for a named detector model, sampled at spacing `deltaF`. Bins at DC and below `fLow` are set to `Infinity`.

**Available models**

*Analytical (LALSim / GWINC formulas)*

| Name | Description |
|---|---|
| `"aLIGOZeroDetHighPower"` | aLIGO zero-detuned, 125 W (standard design sensitivity) |
| `"aLIGOZeroDetLowPower"` | aLIGO zero-detuned, 25 W |
| `"aLIGONSNSOpt"` | aLIGO NS-NS optimized (11° detuning) |
| `"aLIGOBHBH20Deg"` | aLIGO BH-BH 20° detuned, 20 W |
| `"aLIGOHighFrequency"` | aLIGO high-frequency (low SRM transmission) |
| `"AdvVirgo"` | Advanced Virgo design (phenomenological fit) |
| `"KAGRA"` | KAGRA design (phenomenological fit) |
| `"Virgo"` | Initial Virgo |
| `"iLIGOSRD"` | Initial LIGO SRD |

*Tabulated — LIGO-T0900288 / LIGO-P1600143*

| Name | Description |
|---|---|
| `"aLIGOZeroDetHighPowerGWINC"` | aLIGO zero-det high power, GWINC tabulated (LIGO-T0900288) |
| `"EinsteinTelescopeP1600143"` | Einstein Telescope ET-D (LIGO-P1600143) |
| `"CosmicExplorerP1600143"` | Cosmic Explorer baseline (LIGO-P1600143) |
| `"CosmicExplorerPessimisticP1600143"` | CE pessimistic (LIGO-P1600143) |
| `"CosmicExplorerWidebandP1600143"` | CE wideband (LIGO-P1600143) |

*Tabulated — LIGO-T1500293-v13 (Evans, Sturani, Vitale, Hall 2020)*

Design / projected curves:

| Name | Description |
|---|---|
| `"aLIGOT1500293"` | aLIGO broadband (early GWINC run) |
| `"aLIGODesignT1500293"` | aLIGO design sensitivity |
| `"aLIGOPlusT1500293"` | A+ (frequency-dependent squeezing) |
| `"aLIGOPlusSqzOnlyT1500293"` | A+ squeezing upgrade only |
| `"KAGRAT1500293"` | KAGRA design |
| `"KAGRAWidebandT1500293"` | KAGRA wideband configuration |
| `"KAGRASqueezingT1500293"` | KAGRA with squeezing |
| `"EinsteinTelescopeT1500293"` | Einstein Telescope ET-D |
| `"CosmicExplorer1T1500293"` | CE1 (arXiv:1903.04615) |
| `"CosmicExplorer2T1500293"` | CE2 (arXiv:1903.04615) |
| `"AdvancedVirgoT1500293"` | Advanced Virgo design |
| `"AdvancedVirgoWidebandT1500293"` | Advanced Virgo wideband |
| `"AdvancedVirgoSqueezingT1500293"` | Advanced Virgo with squeezing |
| `"VoyagerT1500293"` | Voyager (BlueBird5 configuration) |

Measured / historical curves:

| Name | Description |
|---|---|
| `"LIGOS6T1500293"` | LIGO S6 measured ASD |
| `"LIGOER8T1500293"` | LIGO ER8 (Engineering Run 8) measured ASD |
| `"LIGOO1T1500293"` | LIGO O1 measured ASD |
| `"LIGOO2T1500293"` | LIGO O2 measured ASD |
| `"LIGOH1O3T1500293"` | LIGO H1 O3 measured ASD |
| `"LIGOL1O3T1500293"` | LIGO L1 O3 measured ASD |
| `"VirgoO3T1500293"` | Virgo O3 measured ASD |

*Space-based — LISA (Robson et al. 2019)*

| Name | Description |
|---|---|
| `"LISA"` | LISA (Robson et al. 2019), no confusion noise |
| `"LISAConfusion05yr"` | LISA + galactic confusion noise, 0.5 yr mission |
| `"LISAConfusion1yr"` | LISA + galactic confusion noise, 1 yr |
| `"LISAConfusion2yr"` | LISA + galactic confusion noise, 2 yr |
| `"LISAConfusion4yr"` | LISA + galactic confusion noise, 4 yr |

Returns `$Failed` for unknown names. Use `ListDetectorPSDs[]` to enumerate all names programmatically.

---

### `ListDetectorPSDs`

```mathematica
ListDetectorPSDs[]
```

Returns a list of all available detector PSD model names.

---

### `LoadPSD`

```mathematica
LoadPSD[filename, deltaF, fLow, nFreq]
```

Loads a two-column (frequency [Hz], ASD [1/√Hz]) text file and interpolates onto a regular frequency grid using log-log linear interpolation.

**Arguments**

| Name | Type | Description |
|---|---|---|
| `filename` | string | Path to the text file |
| `deltaF` | real | Target frequency bin spacing [Hz] |
| `fLow` | real | Lower cutoff; bins below this are set to `Infinity` |
| `nFreq` | integer | Number of output frequency bins |

**Returns** a real list of length `nFreq`. Bins beyond the file's frequency range are set to `Infinity`.

---

### `FlatPSD`

```mathematica
FlatPSD[nFreq, deltaF, fLow]
```

Returns a flat (unity) PSD of length `nFreq`. Useful for whitened data or for testing without detector noise. Bins at DC and below `fLow` are set to `Infinity`.

---

## Utilities

### `GetCutoffIndices`

```mathematica
GetCutoffIndices[fLow, fHigh, deltaF, nTime]
```

Returns `{kmin, kmax}` — the 1-based frequency bin indices corresponding to `[fLow, fHigh]`. The convention matches pyCBC's `get_cutoff_indices` (converted from 0-based to 1-based). The summation range in inner products is `kmin` to `kmax - 1` inclusive.

---

### `CyclicTimeShift`

```mathematica
CyclicTimeShift[hTilde, dt, deltaF]
```

Applies a cyclic time shift of `dt` seconds to a frequency-domain waveform by multiplying each bin by $e^{2\pi i f_k \Delta t}$. Does not change the signal norm. Used internally by `GWOptimizedMatch`.

---

## FFT Conventions

GWMMatch uses the **pyCBC/FFTW convention**:

$$\tilde{h}[k] = \Delta t \sum_{j=0}^{N-1} h[j]\, e^{-2\pi i j k / N}$$

In Mathematica this corresponds to `Fourier[h, FourierParameters -> {1, -1}]` scaled by `deltaT`, keeping only the first `Floor[N/2] + 1` bins (one-sided spectrum).

The inverse (used in `GWMatch`) is an unnormalized c2c IDFT:

$$h[j] = \sum_{k=0}^{N-1} \tilde{H}[k]\, e^{+2\pi i j k / N}$$

corresponding to `InverseFourier[H, FourierParameters -> {-1, -1}]`.
