"""Join peak spectra, independent NIST results, and RI-window candidates."""

from bisect import bisect_left, bisect_right
import hashlib
import json
import math
import re

import numpy as np

from .spectra import MAX_MASS, MIN_MASS, cosine, eic_similarity, interpolate, relative


CAS = re.compile(r"[1-9][0-9]{1,6}-[0-9]{2}-[0-9]")
DIAGNOSTIC_IONS = {"112-72-1": (83, 196), "36653-82-4": (83, 224), "112-92-5": (83, 252)}


def sha256(path):
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def local_reference_json(path):
    with path.open(encoding="utf-8-sig") as handle:
        records = json.load(handle)
    if not isinstance(records, list) or not records:
        raise ValueError("Reference JSON must be a nonempty array.")
    for row in records:
        ai = float(row["ai"])
        if not math.isfinite(ai) or ai <= 0 or not row.get("name"):
            raise ValueError("Each reference requires a name and positive AI.")
        row["ai"] = ai
        row.setdefault("cas", "")
        row.setdefault("pdf_page", None)
        row.setdefault("plot_type", "JSON")
        masses = {int(m): float(v) for m, v in row.get("ions", {}).items()}
        if any(not (1 <= m <= 1200 and math.isfinite(v) and v >= 0) for m, v in masses.items()):
            raise ValueError(f"Invalid ions for reference {row['name']}")
        if max(masses.values(), default=0) <= 0:
            row.update(ions={}, spectrum_status="No comparable reference spectrum")
            continue
        lo = int(row.get("mass_min", min(masses)))
        hi = int(row.get("mass_max", max(masses)))
        if not (1 <= lo <= min(masses) <= max(masses) <= hi <= 1200):
            raise ValueError(f"Invalid reference scan limits for {row['name']}")
        row.update(ions=masses, mass_min=lo, mass_max=hi, spectrum_status="Compared")
    return records


def window_matches(peaks, references, queries, links, hits, anchors, window, scan_data,
                   measured_range=(MIN_MASS, MAX_MASS)):
    references = sorted(references, key=lambda entry: entry["ai"])
    indices = [entry["ai"] for entry in references]
    all_matches, all_fragments, diagnostic = [], [], []
    for peak in peaks:
        peak_id = peak["peak_id"]
        modes = links[peak_id]
        ri = interpolate(peak["rt_min"], anchors)
        peak["ri"] = ri
        peak["primary_query"] = modes.get("BG") or modes["RAW"]
        if ri is None:
            continue
        selection = references[bisect_left(indices, ri - window):bisect_right(indices, ri + window)]
        selection = sorted(selection, key=lambda entry: abs(entry["ai"] - ri))
        spectra = {mode: relative(queries[name]) for mode, name in modes.items()}
        mode = "BG" if "BG" in modes else "RAW"
        for entry in selection:
            match_id = len(all_matches) + 1
            cas = entry["cas"]
            same_cas = [f"{k} #{hit['rank']}" for k, query in modes.items() for hit in hits.get(query, [])
                        if CAS.fullmatch(cas) and hit["cas"] == cas]
            lo, hi, score, ion_rows = compare_fragments(entry, spectra, mode, measured_range)
            all_matches.append({
                "reference_id": match_id, "peak_id": peak_id, "sample": peak["sample"],
                "rt_min": peak["rt_min"], "ri": ri, "name": entry["name"], "cas": cas,
                "ai": entry["ai"], "delta_ri": ri - entry["ai"], "pdf_page": entry["pdf_page"],
                "plot_type": entry.get("plot_type", ""), "spectrum_status": entry["spectrum_status"],
                "mass_min": lo, "mass_max": hi, "spectrum": mode, "cosine": score,
                "nist_same_cas": ", ".join(same_cas),
            })
            all_fragments.extend((match_id, peak_id, *row) for row in ion_rows)
            if cas in DIAGNOSTIC_IONS and peak["sample"] in scan_data:
                scans, times = scan_data[peak["sample"]]
                ion, target = DIAGNOSTIC_IONS[cas]
                check = eic_similarity(ion, target, scans, times, peak["rt_min"])
                if check:
                    diagnostic.append((match_id, peak_id, ion, target, *check))
    return all_matches, all_fragments, diagnostic


def compare_fragments(reference, measured, base_mode, measured_range):
    if not reference["ions"]:
        return None, None, None, []
    lower = max(measured_range[0], reference["mass_min"])
    upper = min(measured_range[1], reference["mass_max"])
    if lower > upper:
        return None, None, None, []
    ions = reference["ions"]
    expected = np.array([ions.get(m, 0) for m in range(lower, upper + 1)], dtype=float)
    signal = measured[base_mode][lower:upper + 1]
    similarity = cosine(expected, signal)
    strong = sorted((m for m in ions if lower <= m <= upper), key=lambda m: (-ions[m], m))[:8]
    sample_top = sorted(np.flatnonzero(signal), key=lambda m: (-signal[m], m))[:8]
    keys = set(strong).union(lower + m for m in sample_top)
    if reference["cas"] in DIAGNOSTIC_IONS:
        keys.update(m for m in DIAGNOSTIC_IONS[reference["cas"]] if lower <= m <= upper)
    scale = max(ions.values())
    rows = [(int(m), 100 * ions[m] / scale if m in ions else None,
             float(measured["BG"][m]) if "BG" in measured else None,
             float(measured["RAW"][m]) if "RAW" in measured else None) for m in sorted(keys)]
    return lower, upper, similarity, rows
