"""Peak detection, retention-index interpolation, MSP, and replicate matching."""

from collections import defaultdict
import csv
import json
import math
import re

import numpy as np


PEAK_COLUMNS = ("peak_id", "sample", "rt_min", "tic")
SAFE_ID = re.compile(r"[A-Za-z0-9_.-]+")
PEAK_PAIR = re.compile(r"(\d+(?:\.\d+)?)\s+(\d+(?:\.\d+)?(?:[Ee][-+]?\d+)?)")
MIN_MASS, MAX_MASS = 38, 600
MASS_COUNT = 1201


def read_csv(path, columns):
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        missing = set(columns) - set(reader.fieldnames or ())
        if missing:
            raise ValueError(f"{path.name}: missing fields {', '.join(sorted(missing))}")
        return list(reader)


def write_peak_csv(path, peaks):
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        table = csv.DictWriter(handle, PEAK_COLUMNS, extrasaction="ignore")
        table.writeheader()
        table.writerows(peaks)


def load_anchors(path):
    rows = read_csv(path, ("carbon", "rt_min"))
    anchors = [(int(row["carbon"]), float(row["rt_min"])) for row in rows]
    if len(anchors) < 2:
        raise ValueError("Alkane calibration needs at least two confirmed anchors.")
    for i, (n, rt) in enumerate(anchors):
        if n <= 0 or not math.isfinite(rt) or rt <= 0:
            raise ValueError("Alkane carbon numbers and retention times must be positive.")
        if i and (n != anchors[i - 1][0] + 1 or rt <= anchors[i - 1][1]):
            raise ValueError("Calibration must have consecutive carbon numbers and increasing retention times.")
    return anchors


def interpolate(rt, anchors):
    carbon, times = zip(*anchors)
    if times[0] <= rt <= times[-1]:
        return float(np.interp(rt, times, np.asarray(carbon) * 100))
    return None


def peak_table(path):
    peaks = []
    seen = set()
    for row in read_csv(path, ("peak_id", "sample", "rt_min")):
        name = row["peak_id"]
        if not name or not SAFE_ID.fullmatch(name) or name in seen:
            raise ValueError(f"Peak ID must be unique and ASCII letters/digits/period/underscore/hyphen: {name!r}")
        if not row["sample"]:
            raise ValueError(f"Sample label missing for {name}.")
        row["rt_min"] = float(row["rt_min"])
        if row["rt_min"] <= 0 or not math.isfinite(row["rt_min"]):
            raise ValueError(f"Invalid RT for {name}.")
        row["tic"] = float(row["tic"]) if row.get("tic") else None
        seen.add(name)
        peaks.append(row)
    if not peaks:
        raise ValueError("Peak table is empty.")
    return peaks


def parse_msp(path):
    text = path.read_text(encoding="utf-8-sig")
    result = {}
    for block in re.split(r"(?im)^Name\s*:\s*", text)[1:]:
        name, _, remainder = block.partition("\n")
        name = name.strip()
        if not SAFE_ID.fullmatch(name) or name in result:
            raise ValueError(f"Invalid or duplicate MSP name: {name!r}")
        header = re.search(r"(?im)^Num\s*Peaks\s*:\s*(\d+)\s*$", remainder)
        if not header:
            raise ValueError(f"Missing Num Peaks in MSP record {name}")
        data = remainder[header.end():]
        pairs = PEAK_PAIR.findall(data)
        if len(pairs) != int(header[1]) or not pairs or re.sub(r"[\s;]+", "", PEAK_PAIR.sub("", data)):
            raise ValueError(f"Incorrect or empty Num Peaks in MSP record {name}")
        spectrum = np.zeros(MASS_COUNT)
        for mz, abundance in pairs:
            nominal = math.floor(float(mz) + 0.5)
            amount = float(abundance)
            if not (0 <= nominal < MASS_COUNT and math.isfinite(amount)):
                raise ValueError(f"Invalid mass or abundance in MSP record {name}")
            spectrum[nominal] += amount
        if not spectrum.any():
            raise ValueError(f"MSP record has no positive intensities: {name}")
        result[name] = spectrum
    if not result:
        raise ValueError("The MSP file contains no spectra.")
    return result


def link_spectra(peaks, queries):
    names = {}
    expected = set()
    for peak in peaks:
        name = peak["peak_id"]
        raw, bg = name + "_RAW", name + "_BG"
        if name in queries:
            if raw in queries or bg in queries:
                raise ValueError(f"Use either {name} or {name}_RAW/{name}_BG, not both.")
            names[name] = {"RAW": name}
        else:
            names[name] = {mode: query for mode, query in (("RAW", raw), ("BG", bg)) if query in queries}
        if not names[name]:
            raise ValueError(f"No MSP spectrum for peak {name}")
        expected.update(names[name].values())
    if expected != queries.keys():
        raise ValueError(f"MSP entries without matching peaks: {sorted(queries.keys() - expected)[:8]}")
    return names


def write_msp(path, queries):
    with path.open("w", encoding="ascii", newline="\n") as handle:
        for name, values in queries.items():
            masses = np.flatnonzero(values)
            handle.write(f"Name: {name}\nNum Peaks: {len(masses)}\n")
            for ion in masses:
                handle.write(f"{ion} {values[ion]:.4f}\n")
            handle.write("\n")


def read_scans(path):
    with path.open(encoding="utf-8") as handle:
        scans = [json.loads(line) for line in handle]
    if len(scans) < 15:
        raise ValueError(f"Insufficient scans: {path.name}")
    times = np.array([row["rt"] for row in scans], dtype=float)
    intensities = np.array([sum(row["intensity"]) for row in scans], dtype=float)
    if not np.all(np.isfinite(times)) or not np.all(np.diff(times) > 0):
        raise ValueError(f"Scan RT must increase monotonically: {path.name}")
    return scans, times, intensities


def mean_spectrum(scans, index):
    values = np.zeros(MASS_COUNT)
    for record in scans[index - 1:index + 2]:
        mz = np.floor(np.array(record["mz"]) + 0.5).astype(int)
        intensity = np.asarray(record["intensity"], dtype=float)
        keep = (mz >= 0) & (mz < MASS_COUNT)
        np.add.at(values, mz[keep], intensity[keep] / 3)
    return values


def detect(path, sample_id, sample_name):
    scans, times, tic = read_scans(path)
    smooth = np.convolve(tic, np.ones(5) / 5, "same")
    sampling = np.median(np.diff(times))
    window, gap = int(0.12 / sampling), max(2, int(0.03 / sampling))
    if window <= gap + 3 or len(times) < 2 * window + 8:
        raise ValueError("The scan density is insufficient for the configured TIC and background windows.")
    floor = 0.001 * max(smooth)
    candidates = []
    for i in range(window + 2, len(times) - window - 2):
        if not (smooth[i - 1] < smooth[i] >= smooth[i + 1]):
            continue
        prominence = smooth[i] - max(min(smooth[i - window:i]), min(smooth[i + 1:i + window + 1]))
        change = np.diff(tic[i - window:i + window + 1])
        noise = 1.4826 * np.median(np.abs(change - np.median(change)))
        if prominence > max(floor, 8 * noise):
            candidates.append((prominence, i))
    selected = []
    for _, i in sorted(candidates, reverse=True):
        if all(abs(times[i] - times[j]) >= 0.06 for j in selected):
            selected.append(i)
    peaks, queries = [], {}
    for number, i in enumerate(sorted(selected), start=1):
        i = max(range(i - 2, i + 3), key=lambda j: tic[j])
        name = f"{sample_id}_P{number:03d}"
        left = min(range(i - window, i - gap), key=lambda j: smooth[j])
        right = min(range(i + gap, i + window), key=lambda j: smooth[j])
        raw = mean_spectrum(scans, i)
        bg = np.maximum(raw - (mean_spectrum(scans, left) + mean_spectrum(scans, right)) / 2, 0)
        if not raw.any():
            continue
        peaks.append(dict(peak_id=name, sample=sample_name, rt_min=float(times[i]), tic=float(tic[i])))
        queries[name + "_RAW"] = raw
        if bg.any():
            queries[name + "_BG"] = bg
    return peaks, queries


def cosine(a, b):
    norm = np.linalg.norm(a) * np.linalg.norm(b)
    return float(np.clip(a @ b / norm, 0, 1)) if norm else None


def relative(values):
    scale = max(values)
    return 100 * values / scale if scale else values


def replicate_pairs(peaks, queries, links):
    by_sample = defaultdict(list)
    for peak in peaks:
        by_sample[peak["sample"]].append(peak)
    names = sorted(by_sample)
    results = []
    for i, a in enumerate(names):
        for b in names[i + 1:]:
            possible = []
            for p in by_sample[a]:
                id1 = p["peak_id"]
                p_link = links[id1]
                v1 = queries[p_link.get("BG") or p_link["RAW"]]
                for q in by_sample[b]:
                    id2 = q["peak_id"]
                    delta = q["rt_min"] - p["rt_min"]
                    if abs(delta) > 0.10:
                        continue
                    q_link = links[id2]
                    v2 = queries[q_link.get("BG") or q_link["RAW"]]
                    similarity = cosine(v1, v2)
                    if similarity is not None and similarity >= 0.80:
                        possible.append((similarity - abs(delta), id1, id2, delta, similarity))
            used = set()
            for _, id1, id2, delta, similarity in sorted(possible, reverse=True):
                if id1 in used or id2 in used:
                    continue
                used.update((id1, id2))
                results.append((a, b, id1, id2, delta, similarity))
    return results


def eic_similarity(ion, target, scans, times, rt):
    window = np.flatnonzero(np.abs(times - rt) <= 0.07)
    if len(window) < 7:
        return None
    traces = [[], []]
    for i in window:
        mass = np.floor(np.array(scans[i]["mz"]) + 0.5).astype(int)
        values = np.asarray(scans[i]["intensity"], dtype=float)
        for trace, mz in zip(traces, (ion, target), strict=True):
            trace.append(float(values[mass == mz].sum()))
    smooth = [np.convolve(trace, np.ones(5) / 5, "valid") for trace in traces]
    if any(np.std(trace) == 0 for trace in smooth):
        return None
    local_times = times[window][2:-2]
    focus = abs(local_times - rt) <= 0.03
    if not np.any(focus):
        return None
    apices = [local_times[focus][np.argmax(trace[focus])] for trace in smooth]
    return float(np.corrcoef(*smooth)[0, 1]), float((apices[1] - apices[0]) * 60)
