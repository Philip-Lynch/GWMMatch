# GWMMatch

A gravitational wave match/mismatch package for Wolfram Mathematica.  
GWMMatch computes the noise-weighted match between two waveforms, maximised over time and phase shifts, using a port of the [pyCBC](https://pycbc.org) `matched_filter` / `optimized_match` pipeline.

**Features**

- Noise-weighted inner product, overlap, match, and mismatch
- Discrete (IFFT-based) and Brent-optimised sub-sample match
- Planck taper (port of LALSuite `XLALSimInspiralREAL8WaveTaper`)
- 20+ built-in detector PSD models:
  - Analytical: aLIGO (5 configurations), Advanced Virgo, KAGRA, initial LIGO/Virgo
  - Tabulated: ET-D, CE baseline / pessimistic / wideband (LIGO-P1600143), aLIGO GWINC
  - Space-based: LISA with/without galactic confusion noise (Robson et al. 2019)
- Generic `LoadPSD` for arbitrary two-column (frequency, ASD) files

---

## Installation

Copy the `GWMMatch/` folder into your Mathematica user application directory:

```mathematica
(* Find the right path *)
FileNameJoin[{$UserBaseDirectory, "Applications"}]
```

Then load the package in any notebook or script:

```mathematica
<< GWMMatch`
```

Alternatively, load it directly with an explicit path:

```mathematica
Get["/path/to/GWMMatch/GWMMatch.wl"]
```

---

## Quick Start

```mathematica
<< GWMMatch`

(* Parameters *)
sampleRate = 2048.0;  duration = 2.0;
nTime = Round[sampleRate * duration];
deltaT = 1.0 / sampleRate;  deltaF = 1.0 / duration;
fLow = 20.0;  fHigh = 512.0;
nFreq = Floor[nTime / 2] + 1;

(* Build two waveforms *)
tVals = Table[n deltaT, {n, 0, nTime - 1}];
h1 = GWTaper[Table[Sin[2 Pi 60.0 t], {t, tVals}]];
h2 = GWTaper[RotateLeft[h1, 40]];  (* 40-sample time shift *)

(* Load a detector PSD *)
psd = DetectorPSD["aLIGOZeroDetHighPower", nFreq, deltaF, fLow];

(* Compute the match *)
{match, timeIndex} = GWOptimizedMatch[h1, h2, deltaT, psd, deltaF, fLow, fHigh];
Print["Match = ", N[match, 8]]
Print["Mismatch = ", N[1 - match, 8]]
```

---

## Public API

| Function | Description |
|---|---|
| `GWInnerProduct[h1T, h2T, psd, dF, fL, fH]` | Noise-weighted inner product `<h1\|h2>` |
| `GWSigmaSq[hT, psd, dF, fL, fH]` | Signal norm squared `<h\|h>` |
| `GWSigma[hT, psd, dF, fL, fH]` | Signal norm `Sqrt[<h\|h>]` |
| `GWOverlap[h1T, h2T, psd, dF, fL, fH]` | Normalized overlap (no maximization) |
| `GWMatch[h1, h2, dT, psd, dF, fL, fH]` | Match via IFFT, returns `{match, peakIndex}` |
| `GWOptimizedMatch[h1, h2, dT, psd, dF, fL, fH]` | Sub-sample match via Brent's method |
| `GWMismatch[h1, h2, dT, psd, dF, fL, fH]` | `1 - GWOptimizedMatch`, returns `{mismatch, peakIndex}` |
| `GWTaper[h, method, samples]` | Planck taper; method ∈ `"TAPER_START"`, `"TAPER_END"`, `"TAPER_STARTEND"` |
| `DetectorPSD[name, nFreq, dF, fLow]` | Sample a named PSD onto a regular grid |
| `ListDetectorPSDs[]` | List all available PSD model names |
| `LoadPSD[file, dF, fLow, nFreq]` | Load a two-column (freq, ASD) file |
| `FlatPSD[nFreq, dF, fLow]` | Flat (unity) PSD |
| `GetCutoffIndices[fLow, fHigh, dF, nTime]` | Frequency band → 1-based array indices |
| `CyclicTimeShift[hT, dt, dF]` | Multiply `h(f)` by `exp(2πi f dt)` |

See [docs/API.md](docs/API.md) for full parameter descriptions and examples.

---

## Examples

The `notebooks/` directory contains runnable Wolfram Script (`.wls`) files:

| File | Description |
|---|---|
| [`01-BasicUsage.wls`](notebooks/01-BasicUsage.wls) | End-to-end match between two chirp waveforms |
| [`02-PSDs.wls`](notebooks/02-PSDs.wls) | Plotting and comparing all built-in PSD models |
| [`03-Tapering.wls`](notebooks/03-Tapering.wls) | Planck taper: time-domain and spectral leakage effects |

Run any example from the command line:

```bash
wolframscript -file notebooks/01-BasicUsage.wls
```

---

## Tests

Unit tests live in `tests/` as `.wlt` files and use Mathematica's built-in `TestReport` framework.

```mathematica
(* Run all tests from a notebook *)
<< GWMMatch`
TestReport["tests/TestInnerProduct.wlt"]
TestReport["tests/TestMatch.wlt"]
TestReport["tests/TestPSD.wlt"]
```

Or run everything from the command line:

```bash
wolframscript -file tests/RunTests.wls
```

---

## References

- **pyCBC**: Nitz et al., [pycbc.filter.matchedfilter](https://pycbc.org)
- **LALSuite taper**: McKechan, Robinson & Sathyaprakash, CQG 27 (2010) 084020
- **aLIGO quantum noise**: Buonanno & Chen, PRD 64, 042006 (2001)
- **Next-generation PSDs**: Evans et al., LIGO-P1600143
- **LISA sensitivity**: Robson, Cornish & Liu, CQG 36 (2019) 105011, [arXiv:1803.01944](https://arxiv.org/abs/1803.01944)

---

## License

[MIT](LICENSE)

