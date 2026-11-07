# Peak detection and tracking in noisy XRD scans

This interview-preparation project implements a classical pipeline for locating diffraction peaks and assigning consistent IDs across successive scans. It uses synthetic diffraction-like data with known peak centers and identities.

## Run

From this repository, using Python 3.10 or later:

```sh
python -m pip install -r requirements.txt
python -m unittest -v
python demo.py
```

The command creates `results/metrics.json` and, for each scenario, an example plot, an NPZ of simulated scans, and CSVs of ground-truth peaks and estimated tracks. `benchmark/` contains a checked-in run. To process an external scan, pass increasing, uniformly sampled 2θ angles and finite intensities to `detect(angles, intensities)`. Call `Tracker.update(frame, centers)` once per scan in frame order. The constants in `detect` are calibrated to the demo's angular sampling, peak widths and intensity scale.

## How it works

1. **Smooth and estimate background.** An 11-sample, degree-three Savitzky–Golay filter fits a local polynomial. Grayscale opening takes a running minimum followed by a running maximum over a window about 1.1° wide. This suppresses positive features narrower than the window to estimate a slowly varying background. Broad peaks and scan boundaries can bias this estimate.
2. **Find candidates and refine centers.** SciPy finds local maxima in the smoothed, background-subtracted signal and filters by prominence and spacing. Within ±0.19° of each candidate, bounded least squares fits the unsmoothed, background-subtracted intensities. The pseudo-Voigt is a weighted sum of Gaussian and Lorentzian profiles that share a center and full width at half maximum. The local fit estimates four parameters: height, center, width, and mixture fraction. A fit that reaches its evaluation limit falls back to the candidate location.
3. **Associate across scans.** The assignment matrix uses absolute angular distances. A high cost for pairs beyond 0.12° makes the solver maximize feasible matches before minimizing total distance. Unmatched detections start IDs, and an ID can bridge up to two missing frames. SciPy solves the same linear assignment problem commonly taught with the Hungarian algorithm using its modified Jonker–Volgenant implementation.

## Benchmark protocol

Seeds 21–25 generate 80 scans each per scenario, with 651 points per scan, a curved background, additive Gaussian noise, asymmetric peaks, drift and peak appearance/disappearance. The overlap case adds a crossing peak and higher noise. The same thresholds are used in both cases. The generator and detector were developed together. These runs characterize this simulation and do not estimate field accuracy.

Precision and recall use one-to-one matching within 0.08° per scan. Center MAE uses only matched peaks. Dominant-ID coverage counts, for each true peak in each sequence, how many observations matched its most frequently assigned ID, divided by all true observations. This penalizes gaps and fragmentation and can miss identity merges. Runtime includes detection and tracking and excludes generation, scoring, and plotting.

A single-peak fit cannot reliably separate a merged maximum. A position-only tracker can swap IDs when tracks cross. Phase identification would require interpreting multiple reflections and comparing them with reference patterns under the relevant instrument and material conditions.

[Peak finding](https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.find_peaks.html) · [Linear assignment](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.linear_sum_assignment.html)

