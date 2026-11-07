"""Fit and track peaks in uniformly sampled diffraction scans."""
import numpy as np
from scipy.ndimage import grey_opening
from scipy.optimize import curve_fit, linear_sum_assignment
from scipy.signal import find_peaks, savgol_filter


def profile(x, height, center, width, eta):
    """Height-normalized pseudo-Voigt with shared FWHM."""
    z = 2 * (x - center) / width
    return height * ((1 - eta) * np.exp(-np.log(2) * z**2) + eta / (1 + z**2))


def detect(x, y):
    """Return fitted centers and the estimated baseline. Angles are in degrees."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    if (x.ndim != 1 or x.shape != y.shape or len(x) < 11
            or not np.isfinite(x).all() or not np.isfinite(y).all()):
        raise ValueError("Expected equal finite vectors with at least 11 samples")
    dx = x[1] - x[0]
    if dx <= 0 or not np.allclose(np.diff(x), dx):
        raise ValueError("Expected an increasing uniform angle grid")
    smooth = savgol_filter(y, 11, 3)
    baseline = grey_opening(smooth, size=max(3, 2 * int(0.55 / dx) + 1))
    peaks, _ = find_peaks(smooth - baseline, prominence=9,
                          distance=max(1, round(0.17 / dx)))
    centers = []
    for i in peaks:
        region = np.abs(x - x[i]) <= 0.19
        if region.sum() < 5:
            continue
        try:
            params, _ = curve_fit(
                profile, x[region], (y - baseline)[region],
                p0=[smooth[i] - baseline[i], x[i], 0.15, 0.5],
                bounds=([0, x[i] - 0.08, 0.045, 0],
                        [np.inf, x[i] + 0.08, 0.40, 1]),
                maxfev=300)
            centers.append(params[1])
        except RuntimeError:
            centers.append(x[i])
    return np.array(sorted(centers)), baseline


def match(left, right, gate):
    """Maximize in-gate matches, then minimize total absolute distance."""
    if not len(left) or not len(right):
        return []
    distances = np.abs(np.subtract.outer(left, right))
    penalty = (min(distances.shape) + 1) * gate
    rows, cols = linear_sum_assignment(np.where(distances <= gate, distances, penalty))
    return [(i, j) for i, j in zip(rows, cols) if distances[i, j] <= gate]


class Tracker:
    def __init__(self, gate=0.12, max_gap=2):
        self.gate = gate
        self.max_gap = max_gap
        self.last = {}
        self.next_id = 0

    def update(self, frame, centers):
        """Call once per increasing frame. Return IDs in detection order."""
        self.last = {key: value for key, value in self.last.items()
                     if frame - value[0] <= self.max_gap + 1}
        keys = list(self.last)
        pairs = match([self.last[key][1] for key in keys], centers, self.gate)
        ids = {j: keys[i] for i, j in pairs}
        for j, center in enumerate(centers):
            if j not in ids:
                ids[j] = self.next_id
                self.next_id += 1
            self.last[ids[j]] = (frame, center)
        return [ids[j] for j in range(len(centers))]

