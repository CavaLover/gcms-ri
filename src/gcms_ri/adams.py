"""Read the Adams 2017 ed. 4.1 Appendix III layout from a local PDF."""

from bisect import bisect_left, bisect_right
from io import BytesIO
import math
import re

import numpy as np


HEADER = re.compile(r"\bRT:\s*([\d.]+)\s+AI:\s*(\d+)\s+KI:\s*\d+\s+(.+?)\s+CAS#:\s*([\d-]+)", re.S)

# The raster plots use either 10 or 20 m/z minor ticks. These Adams 4.1
# entries were checked against their printed tick labels and labeled bases.
RASTER_AXES = {
    (1144, "105-21-5"): (10, (85,)), (1159, "7299-40-3"): (10, (43,)),
    (1165, "507-70-0"): (10, (95,)), (1168, "24168-70-5"): (10, (138,)),
    (1187, "21391-84-4"): (10, (41,)), (1191, "41519-18-0"): (10, (55, 70)),
    (1197, "515-88-8"): (10, (85,)), (1219, "77318-48-0"): (10, (43,)),
    (1616, "18312-21-5"): (10, (192,)), (1656, "26184-88-3"): (10, (43,)),
    (1676, "91050-14-5"): (10, (41,)), (1716, "85611-33-2"): (10, (41,)),
    (2387, "90686-31-0"): (20, (81,)), (2387, "70387-35-8"): (20, (81,)),
    (2390, "534-43-0"): (20, (164,)), (2411, "126061-37-8"): (20, (286,)),
    (2422, "10267-31-9"): (20, (43,)), (2443, "3310-97-2"): (20, (135,)),
    (2445, "6811-52-5"): (20, (227,)), (2468, "640-42-6"): (20, (135,)),
    (2500, "629-99-2"): (20, (57,)), (2600, "630-01-3"): (20, (57,)),
    (2700, "593-49-7"): (20, (57,)), (2800, "630-02-4"): (20, (57,)),
    (2900, "630-03-5"): (20, (57,)), (3400, "14167-59-0"): (20, (57,)),
}


def pdf_entries(path, indices, window, index=None):
    try:
        import pdfplumber
    except ImportError as exc:
        raise ValueError("Install the adams extra to read a PDF: pip install 'gcms-ri[adams]'") from exc
    indices = sorted({ri for ri in indices if ri is not None})
    index_count = 0
    entries = []
    with pdfplumber.open(path) as book:
        for page in book.pages:
            try:
                text = page.extract_text(x_tolerance=1.5) or ""
                headers = list(HEADER.finditer(text))
                if not headers:
                    continue
                refs = [
                    {"name": " ".join(hit[3].split()), "cas": hit[4], "ai": int(hit[2]),
                     "pdf_page": page.page_number, "reference_rt_min": float(hit[1])}
                    for hit in headers
                ]
                index_count += len(refs)
                if index is not None:
                    index.extend(dict(r) for r in refs)
                active = [r for r in refs if bisect_left(indices, r["ai"] - window) <
                          bisect_right(indices, r["ai"] + window)]
                if not active:
                    continue
                try:
                    plots = axes(page)
                    if len(plots) != len(refs):
                        raise ValueError(f"{len(plots)} plots for {len(refs)} headers")
                except ValueError as exc:
                    for entry in active:
                        entry.update(ions={}, spectrum_status=f"Review PDF plot: {exc}")
                    entries.extend(active)
                    continue
                words = page.extract_words(x_tolerance=1.1)
                for entry, (_, kind, plot) in zip(refs, plots, strict=True):
                    if entry not in active:
                        continue
                    try:
                        if kind == "vector":
                            ions, lo, hi = vector_spectrum(plot, page.lines, words)
                        else:
                            ions, lo, hi = raster_spectrum(plot, entry)
                        entry.update(ions=ions, mass_min=lo, mass_max=hi,
                                     spectrum_status="Compared", plot_type=kind)
                    except (OSError, IndexError, ValueError) as exc:
                        entry.update(ions={}, spectrum_status=f"Review PDF plot: {exc}", plot_type=kind)
                    entries.append(entry)
            finally:
                page.close()
    if not index_count:
        raise ValueError("No Adams Appendix III records found. Supply the 2017 ed. 4.1 text PDF or --reference-json.")
    return entries, index_count


def near(a, b, eps=0.09):
    return abs(a - b) <= eps


def axes(page):
    out = []
    for curve in page.curves:
        path = curve.get("path", ())
        if [p[0] for p in path] != ["m", "l", "l"]:
            continue
        top, origin, end = (p[1] for p in path)
        if (near(top[0], origin[0]) and near(origin[1], end[1]) and
                end[0] - origin[0] > 250 and origin[1] - top[1] > 50):
            out.append((origin[1], "vector", (top, origin, end)))
    for img in page.images:
        if img["width"] > 250 and 80 < img["height"] < 170:
            out.append(((img["top"] + img["bottom"]) / 2, "raster", img))
    return sorted(out, key=lambda plot: plot[0])


def vector_spectrum(plot, lines, words):
    top, origin, end = plot
    left, base = origin
    right, height = end[0], base - top[1]
    ticks, sticks = [], []
    for line in lines:
        x = line["x0"]
        y0, y1 = line["top"], line["bottom"]
        if not near(x, line["x1"]) or not (left + 0.12 < x < right + 0.09):
            continue
        if near(y0, base) and 1 <= y1 - base <= 6:
            ticks.append(x)
        if near(y1, base) and top[1] - 0.12 <= y0 < base - 0.1:
            sticks.append((x, y0, min(100, 100 * (base - y0) / height)))
    if not ticks:
        raise ValueError("no m/z axis ticks")
    labels = []
    for word in words:
        center = (word["x0"] + word["x1"]) / 2
        if not word["text"].isdigit() or not (base + 0.7 <= word["top"] <= base + 8):
            continue
        x = min(ticks, key=lambda value: abs(value - center))
        if abs(x - center) < 2:
            labels.append((int(word["text"]), x))
    if len(set(labels)) < 4:
        raise ValueError("too few readable m/z tick labels")
    mass, pixels = map(np.array, zip(*sorted(set(labels))))
    slope, offset = np.polyfit(mass, pixels, 1)
    if slope < 0.5 or max(abs(slope * mass + offset - pixels)) > 0.55:
        raise ValueError("inconsistent m/z axis labels")
    annotated = {}
    for word in words:
        if not word["text"].isdigit():
            continue
        label = int(word["text"])
        x = slope * label + offset
        if abs(x - (word["x0"] + word["x1"]) / 2) > 3:
            continue
        options = [stick for stick in sticks if abs(stick[0] - x) < 0.8 * slope
                   and -4 < stick[1] - word["bottom"] < 5]
        if options:
            sx, _, _ = min(options, key=lambda item: abs(item[0] - x))
            annotated.setdefault(sx, set()).add(label)
    ions = {}
    for x, _, abundance in sticks:
        labels_here = annotated.get(x, ())
        m = next(iter(labels_here)) if len(labels_here) == 1 else round((x - offset) / slope)
        ions[m] = max(ions.get(m, 0), abundance)
    if max(ions.values(), default=0) < 98 or len(annotated) < 3:
        raise ValueError("spectrum base peak or annotated ions could not be verified")
    return ions, math.ceil((left - offset) / slope - 0.2), math.floor((right - offset) / slope + 0.2)


def clusters(points):
    if len(points) == 0:
        return []
    return np.split(points, np.where(np.diff(points) != 1)[0] + 1)


def raster_spectrum(img, entry):
    from PIL import Image

    settings = RASTER_AXES.get((entry["ai"], entry["cas"]))
    if settings is None:
        raise ValueError("raster m/z axis not calibrated for this reference")
    minor_step, base_labels = settings
    content = img["stream"].get_data()
    image = np.asarray(Image.open(BytesIO(content)).convert("L"))
    dark = image < 125
    baseline = int(np.argmax(dark.sum(axis=1)))
    row = np.flatnonzero(dark[baseline])
    left, right = int(row[0]), int(row[-1])
    axis_top = baseline - np.flatnonzero(~dark[:baseline, left][::-1])[0]
    lower = dark[baseline + 1:baseline + 7, left:right + 1]
    tick_columns = np.flatnonzero(lower.sum(axis=0) >= 4) + left
    ticks = np.array([np.mean(c) for c in clusters(tick_columns)])
    if len(ticks) < 8 or right - left < 500 or baseline - axis_top < 160:
        raise ValueError("unexpected raster axis")
    masses = 40 + minor_step * np.arange(len(ticks))
    slope, offset = np.polyfit(masses, ticks, 1)
    if max(abs(slope * masses + offset - ticks)) > 1.3:
        raise ValueError("inconsistent raster tick spacing")
    start = left + 2
    height = baseline - axis_top
    patch = dark[axis_top:baseline, start:right + 1][::-1]
    stem_height = np.cumprod(patch[2:], axis=0).sum(axis=0)
    stem_height += 2 * (stem_height > 0)
    ions = {}
    for cols in clusters(np.flatnonzero(stem_height)):
        h = max(stem_height[cols])
        centers = cols[stem_height[cols] >= h - 1]
        mass = (np.mean(centers) + start - offset) / slope
        base = min(base_labels, key=lambda label: abs(label - mass))
        m = base if h >= 0.95 * height and abs(base - mass) < 1 else round(mass)
        ions[m] = max(ions.get(m, 0), min(100, 100.0 * h / height))
    if not any(ions.get(m, 0) >= 96 for m in base_labels):
        raise ValueError("labeled raster base peak could not be verified")
    return ions, math.ceil((left - offset) / slope), math.floor((right - offset) / slope)
