"""Reproduce an experimental opXRD subset and run the existing peak pipeline."""
import argparse
import csv
import hashlib
import json
import shutil
import tempfile
import zipfile
from collections import defaultdict
from pathlib import Path
from time import perf_counter
from urllib.request import urlopen

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from xrd import Tracker, detect

SOURCE_URL = "https://zenodo.org/api/records/14279434/files/opxrd_in_situ.zip/content"
SOURCE_MD5 = "0ab26d52fedde4dae422ddcce8625df9"
PREFIX = "sutter-fella_heymans_0/"
BUNDLE = Path(__file__).parent / "data" / "opxrd-heymans.zip"


def digest(path, algorithm):
    checksum = hashlib.new(algorithm)
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            checksum.update(chunk)
    return checksum.hexdigest()


def scan_key(name):
    """Sort numerically by source group and within-group index."""
    fields = Path(name).stem.split("_")
    return int(fields[-2]), int(fields[-1])


def prepare(archive, destination=BUNDLE):
    """Copy original JSON bytes into a deterministic, attributed subset archive."""
    if digest(archive, "md5") != SOURCE_MD5:
        raise ValueError("Source archive checksum differs from the pinned Zenodo release")
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_suffix(".tmp")
    manifest = []
    try:
        with zipfile.ZipFile(archive) as source, zipfile.ZipFile(temporary, "w") as target:
            names = sorted((n for n in source.namelist()
                            if n.startswith(PREFIX) and n.endswith(".json")), key=scan_key)
            if len(names) != 342:
                raise ValueError("Expected the complete 342-pattern Heymans subset")
            for name in names:
                raw = source.read(name)
                info = zipfile.ZipInfo(name)
                info.compress_type = zipfile.ZIP_DEFLATED
                target.writestr(info, raw)
                manifest.append(dict(member=name, sha256=hashlib.sha256(raw).hexdigest()))
        temporary.replace(destination)
    finally:
        temporary.unlink(missing_ok=True)
    provenance = dict(
        doi="10.5281/zenodo.14279434", source_url=SOURCE_URL, source_md5=SOURCE_MD5,
        license="CC-BY-4.0", subset_sha256=digest(destination, "sha256"),
        selection="All JSON patterns under sutter-fella_heymans_0/", members=manifest)
    destination.with_suffix(".manifest.json").write_text(
        json.dumps(provenance, indent=2) + "\n", encoding="utf-8")


def read_scans(bundle):
    """Validate input and preserve the original uniformly sampled angle grid."""
    manifest = json.loads(bundle.with_suffix(".manifest.json").read_text(encoding="utf-8"))
    if digest(bundle, "sha256") != manifest["subset_sha256"]:
        raise ValueError("Subset checksum mismatch")
    with zipfile.ZipFile(bundle) as archive:
        for name in sorted(archive.namelist(), key=scan_key):
            record = json.loads(archive.read(name))
            x = np.asarray(record["two_theta_values"], float)
            y = np.asarray(record["intensities"], float)
            if (x.ndim != 1 or x.shape != y.shape or len(x) < 11
                    or not np.isfinite(x).all() or not np.isfinite(y).all()
                    or not np.all(np.diff(x) > 0)
                    or not np.allclose(np.diff(x), x[1] - x[0])):
                raise ValueError(f"Invalid or nonuniform angle grid in {name}")
            label = json.loads(record["label"])
            if label["is_simulated"] is not False:
                raise ValueError(f"Expected experimental pattern in {name}")
            yield name, x, y


def analyze(bundle, out):
    out.mkdir(parents=True, exist_ok=True)
    groups = defaultdict(list)
    for name, x, y in read_scans(bundle):
        group, index = scan_key(name)
        groups[group].append((index, name, x, y))
    counts, rows = [], []
    seconds = 0
    scales = {}
    for group, scans in sorted(groups.items()):
        # One scale per source group preserves relative amplitudes across patterns.
        scale = max(float(np.max(y)) for _, _, _, y in scans)
        if scale <= 0:
            raise ValueError(f"Nonpositive intensity scale in group {group}")
        scales[str(group)] = scale
        tracker = Tracker()
        image_rows, example = [], None
        for index, name, x, raw in scans:
            y = raw * (100 / scale)
            start = perf_counter()
            centers, baseline = detect(x, y)
            ids = tracker.update(index, centers)
            seconds += perf_counter() - start
            counts.append(dict(source=name, group=group, source_index=index, peaks=len(centers)))
            rows.extend((group, index, name, track, float(center))
                        for track, center in zip(ids, centers))
            image_rows.append(y)
            if example is None:
                example = x, y, baseline, centers
        if group == 0:
            x, y, baseline, centers = example
            fig, axes = plt.subplots(1, 2, figsize=(12, 4), layout="constrained")
            axes[0].plot(x, y, lw=.7, label="Experimental intensity (scaled)")
            axes[0].plot(x, baseline, lw=.9, label="Estimated background")
            axes[0].scatter(centers, np.interp(centers, x, y), marker="x", label="Detected peaks")
            axes[0].set(xlabel="2θ (degrees)", ylabel="Scaled intensity", title="Source group 0, first pattern")
            axes[0].legend(fontsize=8)
            axes[1].imshow(image_rows, aspect="auto", origin="lower",
                           extent=(x[0], x[-1], -.5, len(scans)-.5), cmap="Greys")
            for track in sorted({r[3] for r in rows if r[0] == 0}):
                observations = [(r[1], r[4]) for r in rows if r[0] == 0 and r[3] == track]
                indices, positions = zip(*observations)
                axes[1].scatter(positions, indices, s=9)
            axes[1].set(xlabel="2θ (degrees)", ylabel="Source pattern index",
                        title="Associations in source index order")
            fig.savefig(out / "experimental.png", dpi=150)
            plt.close(fig)
    with (out / "detections.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(["group", "source_index", "source_member", "track_id", "center_deg"])
        writer.writerows(rows)
    report = dict(
        dataset_doi="10.5281/zenodo.14279434", subset_sha256=digest(bundle, "sha256"),
        experimental_patterns=len(counts), source_groups=len(groups),
        detected_peak_observations=len(rows), pipeline_seconds=seconds,
        intensity_scaling="100 * raw / maximum raw intensity within source group",
        group_intensity_maxima=scales,
        accuracy=None,
        accuracy_note="No hand-labeled peak centers or track IDs supplied for this subset.",
        ordering_note="Numeric source indices within each group; physical timing and temporal continuity unverified.",
        patterns=counts)
    (out / "summary.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({key: value for key, value in report.items()
                      if key not in ("patterns", "group_intensity_maxima")}, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive", type=Path, help="Rebuild subset from the pinned source ZIP")
    parser.add_argument("--download", action="store_true", help="Download 373 MB and rebuild subset")
    parser.add_argument("--out", type=Path, default=Path("results/public"))
    args = parser.parse_args()
    if args.archive and args.download:
        parser.error("Choose either --archive or --download")
    if args.download:
        with tempfile.TemporaryDirectory() as directory:
            archive = Path(directory) / "source.zip"
            with urlopen(SOURCE_URL, timeout=60) as response, archive.open("wb") as stream:
                shutil.copyfileobj(response, stream)
            prepare(archive)
    elif args.archive:
        prepare(args.archive)
    analyze(BUNDLE, args.out)
