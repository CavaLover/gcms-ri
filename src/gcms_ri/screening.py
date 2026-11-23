"""Blank evidence and tentative homolog assignments."""

from dataclasses import dataclass
import re

from .candidates import top_ions
from .chromatography import ScanProfile, shape_similarity
from .spectra import MASS_COUNT, cosine


SILICON_IONS = {73, 147, 207, 221, 267, 281, 341, 355, 401, 429, 503}
ALKANE_IONS = {43,57,71,85,99,113}
N_ALKANE_CAS = dict(enumerate((
    "1120-21-4", "112-40-3", "629-50-5", "629-59-4", "629-62-9",
    "544-76-3", "629-78-7", "593-45-3", "629-92-5", "112-95-8",
    "629-94-7", "629-97-0", "638-67-5", "646-31-1", "629-99-2",
    "630-01-3", "593-49-7", "630-02-4", "630-03-5", "638-68-6",
    "630-04-6", "544-85-4", "630-05-7", "14167-59-0", "630-07-9", "630-06-8",
), 11))


@dataclass(frozen=True)
class Blank:
    label: str
    scope: str
    peaks: tuple
    spectra: dict
    profile: ScanProfile | None


def putative_background(ions, hits):
    if not ions:
        return None
    silicon = set(ions[:6]) & SILICON_IONS
    names = [hit["name"] for hit in hits]
    if (ions[0] in SILICON_IONS and len(silicon) >= 3) or (
        len(silicon) >= 2 and names and all("siloxan" in n.lower() for n in names)
    ):
        return "硅氧烷样背景"
    if ions[0] == 149 and names and all(
        re.search(r"phthalate|phthalic acid|1,2-benzenedicarboxylic", n, re.I) for n in names
    ):
        return "邻苯二甲酸酯类背景候选"
    return None


def compare_blank(rt, spectra, profile, trace, signal, ions, blank, rt_window):
    blank_trace = blank.profile.eic(ions) if blank.profile else None
    candidates = []
    for peak in blank.peaks:
        if abs(peak["rt_min"]-rt) > rt_window:
            continue
        modes = blank.spectra[peak["peak_id"]]
        mode = "BG" if "BG" in spectra and "BG" in modes else "RAW"
        if mode not in spectra or mode not in modes:
            continue
        score = cosine(spectra[mode], modes[mode])
        if score is not None and score >= .80:
            b = blank.profile.nearby_peak(blank_trace, peak["rt_min"], .032) if blank.profile else None
            candidates.append((score, peak["rt_min"], peak["peak_id"], mode, b))
    if not candidates and blank.profile:
        mode = "BG" if "BG" in spectra else "RAW"
        for i in blank.profile.search(blank_trace, rt, rt_window):
            center = blank.profile.times[i]
            b = blank.profile.measure(blank_trace, center)
            if not b or b.signal_to_noise < 5:
                continue
            score = cosine(spectra[mode], blank.profile.spectrum(b, mode))
            if score is not None and score >= .80:
                candidates.append((score, float(center), None, mode, b))
    if not candidates:
        return None
    score, center, pid, mode, b = max(candidates, key=lambda c: c[0])
    shape = shape_similarity(profile, trace, signal, blank.profile, blank_trace, b) if signal and b else None
    blank_rt = b.rt_min if b else center
    aligned = signal is not None and abs(blank_rt-signal.rt_min) <= rt_window
    ratio = signal.area / b.area if signal and b and b.area > 0 and aligned and shape is not None and shape >= .60 else None
    return dict(blank=blank.label, scope=blank.scope, peak_id=pid, rt_min=blank_rt,
                mode=mode, cosine=score, shape_r=shape, area=b.area if b else None, sample_blank_ratio=ratio)


def screen(peak, spectra, hits, profile, blanks, rt_window):
    spectrum = spectra["BG"] if "BG" in spectra else spectra["RAW"]
    ions = [ion["mz"] for ion in top_ions(spectrum, (1,MASS_COUNT-1))]
    family = putative_background(ions, hits)
    diagnostic = (
        [149,223] if family == "邻苯二甲酸酯类背景候选" and spectrum[223] > 0
        else [mz for mz in ions if mz in SILICON_IONS][:3] if family == "硅氧烷样背景"
        else ions[:3]
    )
    rt = peak["rt_min"]
    trace = profile.eic(diagnostic) if profile and diagnostic else None
    signal = profile.nearby_peak(trace, rt,.032) if trace is not None else None
    matches = [match for blank in blanks if diagnostic and
               (match := compare_blank(rt,spectra,profile,trace,signal,diagnostic,blank,rt_window)) is not None]
    dominated = any(m["sample_blank_ratio"] is not None and m["sample_blank_ratio"] <= 1 for m in matches)
    decision = family or ("CK 主导" if dominated else "CK 同峰待分辨" if matches else "保留待鉴定")
    return dict(decision=decision, excluded=bool(family or matches), diagnostic_ions=diagnostic,
                sample_area=signal.area if signal else None, blank_matches=matches)


def alkane_candidate(peak, spectrum, profile, anchors, rt_window, rt_offset):
    ions = {ion["mz"] for ion in top_ions(spectrum, (1,MASS_COUNT-1),6)}
    if len(ions & ALKANE_IONS) < 3 or len(ions & {57,71,85}) < 2:
        return None
    rt = peak["rt_min"]
    choices = [(carbon,standard_rt) for carbon,standard_rt in anchors
               if abs(rt-standard_rt-rt_offset) <= rt_window]
    if len(choices) != 1:
        return None
    carbon, standard_rt = choices[0]
    parent = 14*carbon + 2
    shift, shape = None, None
    if profile and parent < MASS_COUNT:
        fragments = profile.eic([57,71,85])
        molecular = profile.eic([parent])
        fragment_peak = profile.nearby_peak(fragments,rt,.035)
        if fragment_peak:
            molecular_peak = profile.nearby_peak(molecular,fragment_peak.rt_min,.025)
            if molecular_peak:
                shift = molecular_peak.rt_min-fragment_peak.rt_min
                # On the same absolute time axis, do not shift one ion to fit another.
                shape = shape_similarity(profile,fragments,fragment_peak,profile,molecular,fragment_peak)
    return dict(carbon=carbon, standard_rt_min=standard_rt, parent_mz=parent,
                apex_shift_min=shift, shape_r=shape, supported=shape is not None and shape >= .55,
                cas=N_ALKANE_CAS.get(carbon))
