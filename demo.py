"""Generate synthetic scans, benchmark detections and identities, and save plots."""
import argparse
import csv
import json
from collections import Counter
from pathlib import Path
from time import perf_counter

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from xrd import Tracker, detect, match


def synthetic(seed, stress=False):
    rng = np.random.default_rng(seed)
    x = np.linspace(24, 37, 651)
    spectra, truth = [], []
    for frame in range(80):
        f = frame / 79
        locations = {0: 26 + 0.55 * f, 1: 30.4 - 0.38 * f}
        if 15 <= frame < 70:
            locations[2] = 34 + 0.45 * f
        if stress:
            locations[3] = 29.5 + 1.35 * f
        y = 13 + 0.25 * (x - 29)**2 + 2 * np.sin(x + f)
        for key, center in locations.items():
            # Split-width peaks deliberately differ from the symmetric fitted model.
            width = np.where(x < center, 0.067, 0.085) + 0.006 * key
            z = (x - center) / width
            y += (43 - 5 * key) * (0.78 * np.exp(-np.log(2) * z**2) + 0.22 / (1 + z**2))
        spectra.append(y + rng.normal(0, 2.8 if stress else 1.5, x.size))
        truth.append(locations)
    return x, np.array(spectra), truth


def run(seeds, stress, out):
    found = expected = 0
    errors = []
    identities = Counter()
    elapsed = 0
    for seed in seeds:
        x, spectra, truth = synthetic(seed, stress)
        tracker = Tracker()
        points = []
        start = perf_counter()
        for frame, y in enumerate(spectra):
            centers, baseline = detect(x, y)
            ids = tracker.update(frame, centers)
            points.extend((frame, track, center) for track, center in zip(ids, centers))
        elapsed += perf_counter() - start
        found += len(points)
        expected += sum(map(len, truth))
        for frame, known in enumerate(truth):
            observations = [(track, center) for t, track, center in points if t == frame]
            for i, j in match(list(known.values()), [p[1] for p in observations], 0.08):
                track, center = observations[j]
                key = list(known)[i]
                errors.append(center - known[key])
                identities[seed, key, track] += 1
        if seed == seeds[0]:
            np.savez_compressed(out / "example.npz", angles=x, intensity=spectra)
            with (out / "truth.csv").open("w", newline="") as f:
                writer = csv.writer(f)
                writer.writerow(["frame", "true_peak", "angle_deg"])
                writer.writerows((t, key, center) for t, scan in enumerate(truth)
                                 for key, center in scan.items())
            with (out / "tracks.csv").open("w", newline="") as f:
                writer = csv.writer(f)
                writer.writerow(["frame", "track_id", "angle_deg"])
                writer.writerows(points)
            fig, axes = plt.subplots(1, 2, figsize=(10, 4), layout="constrained")
            axes[0].plot(x, spectra[-1], label="Observed", lw=0.9)
            axes[0].plot(x, baseline, label="Estimated background")
            axes[0].set(xlabel="2θ (degrees)", ylabel="Intensity", title="Final scan")
            axes[0].legend()
            axes[1].imshow(spectra, extent=(x[0], x[-1], 79.5, -0.5),
                           aspect="auto", cmap="Greys")
            for track in sorted({row[1] for row in points}):
                rows = [(t, c) for t, k, c in points if k == track]
                t, c = zip(*rows)
                axes[1].scatter(c, t, s=6)
            axes[1].set(xlabel="2θ (degrees)", ylabel="Frame", title="Fitted tracks")
            fig.savefig(out / "example.png", dpi=140)
            plt.close(fig)
    longest = {}
    for (seed, key, track), count in identities.items():
        longest[seed, key] = max(longest.get((seed, key), 0), count)
    return dict(seeds=list(seeds), scans=80 * len(seeds), truth_peaks=expected,
                detections=found, matches=len(errors), precision=len(errors) / max(found, 1),
                recall=len(errors) / expected,
                center_mae_deg=float(np.mean(np.abs(errors))) if errors else None,
                dominant_id_coverage=sum(longest.values()) / expected,
                pipeline_seconds=elapsed)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=Path("results"))
    parser.add_argument("--start-seed", type=int, default=21)
    parser.add_argument("--sequences", type=int, default=5)
    args = parser.parse_args()
    if args.start_seed < 0 or args.sequences < 1:
        parser.error("Expected a nonnegative seed and at least one sequence")
    report = {}
    for name in ("nominal", "overlap"):
        folder = args.out / name
        folder.mkdir(parents=True, exist_ok=True)
        report[name] = run(range(args.start_seed, args.start_seed + args.sequences),
                           name == "overlap", folder)
    output = json.dumps(report, indent=2)
    (args.out / "metrics.json").write_text(output + "\n", encoding="utf-8")
    print(output)

