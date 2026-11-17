"""Command-line orchestration for GC–MS candidate review."""

import argparse
from datetime import datetime, timezone
from importlib.metadata import PackageNotFoundError, version
import json
import math
from pathlib import Path
import platform
import re
import shutil
import subprocess
import sys
import tempfile

from .adams import pdf_entries
from .agilent import compile_reader, native_extract
from .nist import nist_ready, parse_report, search
from .report import book
from .review import local_reference_json, sha256, window_matches
from .spectra import (detect, interpolate, link_spectra, load_anchors, parse_msp,
                      peak_table, read_scans, replicate_pairs, write_msp, write_peak_csv)


def pkg_version():
    try:
        return version("gcms-ri")
    except PackageNotFoundError:
        return "source"


def absolute(text):
    return Path(text).expanduser().resolve()


def positive(text):
    value = float(text)
    if not math.isfinite(value) or value <= 0:
        raise argparse.ArgumentTypeError("enter a finite positive number")
    return value


def parser():
    app = argparse.ArgumentParser(prog="gcms-ri", description="GC–MS library, retention-index, and fragment review.")
    app.add_argument("--version", action="version", version=f"%(prog)s {pkg_version()}")
    sub = app.add_subparsers(dest="command", required=True)
    run = sub.add_parser("run", help="detect peaks from Agilent .D and perform candidate review")
    selection = run.add_mutually_exclusive_group(required=True)
    selection.add_argument("--sample", type=absolute, action="append", help="repeat for each .D or data.ms")
    selection.add_argument("--sample-dir", type=absolute, help="directory whose immediate children are sample .D folders")
    add_agilent(run)
    add_review(run)
    run.add_argument("--skip-nist", action="store_true", help="run RI and fragment review without a NIST report")
    review = sub.add_parser("review", help="review an exported peak table and matching MSP")
    review.add_argument("--peaks", type=absolute, required=True)
    review.add_argument("--spectra", type=absolute, required=True)
    add_review(review)
    standard = sub.add_parser("standard", help="export a standard's candidate peak times and MSP for manual anchor review")
    standard.add_argument("--standard", type=absolute, required=True)
    standard.add_argument("--out", type=absolute, required=True)
    add_agilent(standard)
    nist = sub.add_parser("nist", help="submit an existing ASCII MSP to NIST MS Search")
    nist.add_argument("msp", type=absolute)
    nist.add_argument("--out", type=absolute, required=True)
    add_nist(nist, required=True)
    doctor = sub.add_parser("doctor", help="check optional local Windows components")
    doctor.add_argument("--agilent-dir", type=absolute)
    doctor.add_argument("--agilent-arch", choices=("x64", "x86"), default="x64")
    doctor.add_argument("--nist-dir", type=absolute)
    demo = sub.add_parser("demo", help="run a fully synthetic example without vendor software")
    demo.add_argument("--out", type=absolute, required=True)
    return app


def add_agilent(app):
    app.add_argument("--agilent-dir", type=absolute, required=True, help="directory containing MSDChemDataFile.dll")
    app.add_argument("--agilent-arch", choices=("x64", "x86"), default="x64")
    app.add_argument("--reader-timeout", type=positive, default=180)


def add_nist(app, required=False):
    group = app if required else app.add_mutually_exclusive_group()
    group.add_argument("--nist-dir", type=absolute, required=required, help="NIST MS Search installation directory")
    if not required:
        group.add_argument("--nist-report", type=absolute, help="previous NIST background results.txt")
    app.add_argument("--nist-timeout", type=positive, default=600)
    app.add_argument("--nist-temp", type=absolute, help="existing ASCII-only staging directory")
    if not required:
        app.add_argument("--nist-library-version", help="library/version shown in NIST About or supplied with an imported report")


def add_review(app):
    app.add_argument("--alkanes", type=absolute, required=True, help="verified carbon,rt_min calibration CSV")
    source = app.add_mutually_exclusive_group(required=True)
    source.add_argument("--adams", type=absolute, help="local Adams 2017 ed. 4.1 PDF")
    source.add_argument("--reference-json", type=absolute, help="local AI and EI spectra records")
    app.add_argument("--ri-window", type=positive, default=20)
    app.add_argument("--mass-min", type=int, default=38, help="first acquired nominal m/z used in fragment comparison")
    app.add_argument("--mass-max", type=int, default=600, help="last acquired nominal m/z used in fragment comparison")
    app.add_argument("--out", type=absolute, required=True)
    add_nist(app)


def create_output(path):
    if path.exists():
        raise ValueError(f"Choose a new output directory: {path}")
    path.mkdir(parents=True)


def stage_native(args, paths):
    paths = list(paths)
    if not paths:
        raise ValueError("No sample .D directories were found.")
    sample_names = [p.stem for p in paths]
    if len(set(sample_names)) != len(paths):
        raise ValueError("Use distinct sample directory or file names.")
    work = args.out / "scans"
    work.mkdir()
    peaks, queries, native_scans, manifest = [], {}, {}, []
    labeled = [(f"s{i}", path) for i, path in enumerate(paths, start=1)]
    for label, original, scans_path in native_extract(labeled, args.agilent_dir, work, args.agilent_arch, args.reader_timeout):
        p, q = detect(scans_path, label, original.stem)
        if not p:
            raise ValueError(f"No peaks found in {original.name}. Check the sample or import an external peak table.")
        peaks.extend(p)
        queries.update(q)
        scans, t, _ = read_scans(scans_path)
        native_scans[original.stem] = (scans, t)
        raw_data = original / "data.ms" if original.is_dir() else original
        manifest.append({"label": label, "name": original.stem, "source": str(raw_data), "sha256": sha256(raw_data),
                         "scan_count": len(scans), "peak_count": len(p)})
        print(f"{label} {original.name}: {len(scans)} scans, {len(p)} peaks.", flush=True)
    write_peak_csv(args.out / "peaks.csv", peaks)
    write_msp(args.out / "queries.msp", queries)
    return peaks, queries, native_scans, manifest


def existing_hits(args, msp, queries):
    if args.nist_dir:
        print(f"Submitting {len(queries)} spectra to NIST MS Search.", flush=True)
        report = search(msp, args.nist_dir, args.out / "nist", args.nist_timeout, args.nist_temp)
    elif args.nist_report:
        report = args.out / "nist" / "results.txt"
        report.parent.mkdir()
        shutil.copyfile(args.nist_report, report)
    else:
        return {}, {}
    hits = parse_report(report)
    if hits.keys() != queries.keys():
        missing = sorted(queries.keys() - hits.keys())
        extra = sorted(hits.keys() - queries.keys())
        raise ValueError(f"NIST report/query mismatch. Missing: {missing[:8]}; extra: {extra[:8]}.")
    info = {"nist_report_sha256": sha256(report)}
    settings = report.with_name("settings_before.ini")
    if settings.exists():
        info["nist_settings_sha256"] = sha256(settings)
        info["nist_executable_sha256"] = sha256(nist_ready(args.nist_dir))
        setting = re.search(rb"(?mi)^Use RI in Search\s*=\s*(\d+)", settings.read_bytes())
        if setting:
            info["nist_use_ri_in_search"] = setting[1].decode("ascii")
    if args.nist_library_version:
        info["nist_library_version"] = args.nist_library_version
    return hits, info


def analyze(args, peaks, queries, native_scans, sources, anchors):
    links = link_spectra(peaks, queries)
    msp = args.out / "queries.msp"
    hits, nist_info = existing_hits(args, msp, queries)
    ri = [interpolate(p["rt_min"], anchors) for p in peaks]
    source = args.adams or args.reference_json
    if args.adams:
        print("Reading Adams Appendix III and digitizing spectra within the RI windows.", flush=True)
        references, ref_count = pdf_entries(source, ri, args.ri_window)
    else:
        references = local_reference_json(source)
        ref_count = len(references)
    matches, fragments, diagnostic = window_matches(
        peaks, references, queries, links, hits, anchors, args.ri_window, native_scans,
        (args.mass_min, args.mass_max),
    )
    pairs = replicate_pairs(peaks, queries, links)
    workbook = args.out / "report.xlsx"
    metadata = {
        "gcms-ri version": pkg_version(), "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "python": platform.python_version(), "ri_window": args.ri_window,
        "reference": str(source), "reference_sha256": sha256(source), "reference_index_count": ref_count,
        "calibration_sha256": sha256(args.alkanes), "query_sha256": sha256(msp),
        "mass_range_nominal": f"{args.mass_min}–{args.mass_max}", "sample_sources": sources,
        "peak_detection": "5-scan TIC moving average, 0.12 min prominence window, prominence > max(0.1% global TIC max, 8× scaled MAD of first differences), 0.06 min separation.",
        "peak_spectra": "Three-scan apex mean. BG subtracts means around both local TIC minima about 0.03–0.12 min from the apex. Negative BG intensities are clipped.",
        "replicate_matching": "Pairwise across samples, within 0.10 min and cosine >= 0.80, one-to-one greedy matching.",
        "eic_method": "Pearson r from five-scan-smoothed nominal EICs within RT ±0.07 min. Apex shift uses RT ±0.03 min.",
        **nist_info,
    }
    if not native_scans:
        metadata["peak_detection"] = "Imported peak table and spectra."
        metadata["peak_spectra"] = "Imported MSP; processing was performed by the source application."
    shutil.copyfile(args.alkanes, args.out / "alkanes.csv")
    partial = args.out / "report.partial.xlsx"
    book(partial, peaks, hits, matches, fragments, pairs, diagnostic, anchors, metadata)
    partial.replace(workbook)
    (args.out / "run.json").write_text(json.dumps(metadata, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"{len(peaks)} peaks, {len(matches)} RI candidates, "
          f"{sum(m['cosine'] is not None for m in matches)} compared reference spectra.")
    print(workbook)


def execute(args):
    if args.command == "doctor":
        print(f"gcms-ri {pkg_version()} | Python {platform.python_version()} | {platform.system()}")
        if args.nist_dir:
            print(f"NIST executable: {nist_ready(args.nist_dir)}")
        if args.agilent_dir:
            with tempfile.TemporaryDirectory(prefix="gcms-ri-doctor-") as temp:
                compile_reader(args.agilent_dir, Path(temp), args.agilent_arch)
            print(f"Agilent reader compiled ({args.agilent_arch}): {args.agilent_dir}")
        return
    if args.command == "demo":
        from .demo import demo
        demo(args.out, execute)
        return
    if args.command == "nist":
        parse_msp(args.msp)
        args.out.parent.mkdir(parents=True, exist_ok=True)
        print(search(args.msp, args.nist_dir, args.out, args.nist_timeout, args.nist_temp))
        return
    if args.command == "standard":
        create_output(args.out)
        stage_native(args, [args.standard])
        print("Calibrate by confirming the alkane identities and entering carbon,rt_min in an anchor CSV.")
        return
    if not (1 <= args.mass_min <= args.mass_max <= 1200):
        raise ValueError("Require 1 <= --mass-min <= --mass-max <= 1200.")
    anchors = load_anchors(args.alkanes)
    if args.command == "run":
        if args.skip_nist and (args.nist_dir or args.nist_report):
            raise ValueError("Choose NIST input or --skip-nist.")
        if not (args.skip_nist or args.nist_dir or args.nist_report):
            raise ValueError("Use --nist-dir, --nist-report, or --skip-nist explicitly for native runs.")
        if args.sample_dir and not args.sample_dir.is_dir():
            raise ValueError(f"Sample directory not found: {args.sample_dir}")
        paths = args.sample or sorted(p for p in args.sample_dir.iterdir() if p.is_dir() and p.suffix.lower() == ".d")
        create_output(args.out)
        peaks, queries, native, sources = stage_native(args, paths)
    else:
        peaks, queries = peak_table(args.peaks), parse_msp(args.spectra)
        link_spectra(peaks, queries)
        create_output(args.out)
        write_peak_csv(args.out / "peaks.csv", peaks)
        write_msp(args.out / "queries.msp", queries)
        native = {}
        sources = [{"peak_table_sha256": sha256(args.peaks), "source_msp_sha256": sha256(args.spectra)}]
    analyze(args, peaks, queries, native, sources, anchors)


def main(argv=None):
    app = parser()
    args = app.parse_args(argv)
    try:
        execute(args)
    except (ValueError, OSError, KeyError, TypeError, subprocess.SubprocessError) as exc:
        print(f"gcms-ri: {exc}", file=sys.stderr)
        return 2
    except KeyboardInterrupt:
        print("gcms-ri: interrupted.", file=sys.stderr)
        return 130
    return 0


if __name__ == "__main__":
    sys.exit(main())
