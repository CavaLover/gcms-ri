"""Local extracted-ion chromatograms used for study review."""

from dataclasses import dataclass

import numpy as np

from .spectra import MASS_COUNT, mean_spectrum, read_scans


@dataclass(frozen=True)
class IonPeak:
    rt_min: float
    height: float
    area: float
    noise: float

    @property
    def signal_to_noise(self):
        return self.height / max(1.0, self.noise)


class ScanProfile:
    """Nominal ions from one scan file, with reversible local integrations."""

    def __init__(self, path):
        self.scans, self.times, _ = read_scans(path)
        counts = [len(scan["mz"]) for scan in self.scans]
        self.indices = np.repeat(np.arange(len(self.scans)), counts)
        self.masses = np.floor(np.concatenate([scan["mz"] for scan in self.scans]) + 0.5).astype(int)
        self.current = np.concatenate([scan["intensity"] for scan in self.scans])

    def eic(self, ions):
        selected = np.isin(self.masses, ions)
        values = np.bincount(self.indices[selected], weights=self.current[selected], minlength=len(self.times))
        return np.convolve(values, np.ones(5) / 5, "same")

    def search(self, trace, center, radius):
        window = np.flatnonzero(abs(self.times - center) <= radius)
        maxima = [i for i in window if 0 < i < len(trace)-1 and trace[i-1] <= trace[i] > trace[i+1]]
        return sorted(maxima, key=lambda i: trace[i], reverse=True)

    def measure(self, trace, center):
        t = self.times
        left = (t >= center-.135) & (t <= center-.095)
        right = (t >= center+.095) & (t <= center+.135)
        within = abs(t - center) <= .065
        if min(np.count_nonzero(left), np.count_nonzero(right)) < 3 or np.count_nonzero(within) < 7:
            return None
        l, r = np.median(trace[left]), np.median(trace[right])
        baseline = np.interp(t[within], [np.median(t[left]), np.median(t[right])], [l,r])
        corrected = np.maximum(trace[within] - baseline, 0)
        # Trapezoid integral on native scan times, including irregular grids.
        area = np.dot(np.diff(t[within]), (corrected[:-1] + corrected[1:]) / 2)
        flanks = np.concatenate((trace[left]-l, trace[right]-r))
        noise = 1.4826 * np.median(abs(flanks - np.median(flanks)))
        return IonPeak(float(center), float(np.max(corrected)), float(area), float(noise))

    def nearby_peak(self, trace, rt, tolerance):
        for i in self.search(trace, rt, tolerance):
            peak = self.measure(trace, self.times[i])
            if peak and peak.signal_to_noise >= 5:
                return peak
        return None

    def spectrum(self, peak, mode):
        fit = int(np.argmin(abs(self.times - peak.rt_min)))
        raw = mean_spectrum(self.scans, fit)
        if mode == "RAW":
            return raw
        left = np.flatnonzero((self.times >= peak.rt_min-.135) & (self.times <= peak.rt_min-.095))
        right = np.flatnonzero((self.times >= peak.rt_min+.095) & (self.times <= peak.rt_min+.135))
        baseline = sum(
            (mean_spectrum(self.scans, int(i)) for i in np.r_[left,right] if 0 < i < len(self.times)-1),
            np.zeros(MASS_COUNT),
        )
        count = np.count_nonzero((np.r_[left,right] > 0) & (np.r_[left,right] < len(self.times)-1))
        return np.maximum(raw - baseline / count, 0) if count else raw


def shape_similarity(profile_a, trace_a, peak_a, profile_b, trace_b, peak_b):
    relative_time = np.linspace(-.036, .036, 25)
    a = np.interp(relative_time+peak_a.rt_min, profile_a.times, trace_a)
    b = np.interp(relative_time+peak_b.rt_min, profile_b.times, trace_b)
    if np.std(a) and np.std(b):
        return float(np.clip(np.corrcoef(a,b)[0,1],-1,1))
    return None
