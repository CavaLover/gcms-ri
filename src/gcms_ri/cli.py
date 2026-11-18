"""Command-line orchestration for GC–MS candidate review."""

import argparse
from datetime import datetime, timezone
from importlib.metadata import PackageNotFoundError, version
from io import TextIOWrapper
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
from .candidate_report import write_reports
from .candidates import summarize
from .evidence import load_evidence
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
    app = argparse.ArgumentParser(prog="gcms-ri", description="GC–MS 谱库、保留指数、离子与文献复核。")
    app.add_argument("--version", action="version", version=f"%(prog)s {pkg_version()}")
    sub = app.add_subparsers(dest="command", required=True)
    run = sub.add_parser("run", help="读取Agilent .D、检索谱库并生成中文候选汇总")
    selection = run.add_mutually_exclusive_group(required=True)
    selection.add_argument("--sample", type=absolute, action="append", help="repeat for each .D or data.ms")
    selection.add_argument("--sample-dir", type=absolute, help="directory whose immediate children are sample .D folders")
    run.add_argument("--blank", type=absolute, action="append", default=[], help="CK的.D路径，可重复；不计为生物学重复")
    add_agilent(run)
    add_review(run)
    run.add_argument("--skip-nist", action="store_true", help="run RI and fragment review without a NIST report")
    review = sub.add_parser("review", help="复核已有峰表、MSP与NIST报告")
    review.add_argument("--peaks", type=absolute, required=True)
    review.add_argument("--spectra", type=absolute, required=True)
    add_review(review)
    summary = sub.add_parser("summarize", help="复用旧分析目录，仅重做中文候选与文献汇总")
    summary.add_argument("--run-dir", type=absolute, required=True)
    summary.add_argument("--reference-index", type=absolute, help="完整CAS/AI索引JSON，不能仅提供RI窗口内子集")
    summary.add_argument("--index-kind", choices=("Adams","local"), default="local", help="外部索引来源，默认不宣称是Adams")
    summary.add_argument("--out", type=absolute, required=True)
    summary.add_argument("--ri-window", type=positive, default=20)
    summary.add_argument("--mass-min", type=int, default=38)
    summary.add_argument("--mass-max", type=int, default=600)
    add_summary(summary)
    literature = sub.add_parser("literature", help="按CAS查询公开NIST RI与文献线索（不上传样品）")
    literature.add_argument("--cas", action="append", required=True)
    literature.add_argument("--out", type=absolute, required=True)
    literature.add_argument("--timeout", type=positive, default=20)
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
    source = app.add_mutually_exclusive_group()
    source.add_argument("--adams", type=absolute, help="local Adams 2017 ed. 4.1 PDF")
    source.add_argument("--reference-json", type=absolute, help="local AI and EI spectra records")
    app.add_argument("--ri-window", type=positive, default=20)
    app.add_argument("--mass-min", type=int, default=38, help="first acquired nominal m/z used in fragment comparison")
    app.add_argument("--mass-max", type=int, default=600, help="last acquired nominal m/z used in fragment comparison")
    app.add_argument("--out", type=absolute, required=True)
    add_nist(app)
    add_summary(app)


def score_threshold(text):
    value=int(text)
    if not 0 <= value < 999:
        raise argparse.ArgumentTypeError('匹配分数阈值须在0–998之间。')
    return value


def add_summary(app):
    app.add_argument('--mf-min',type=score_threshold,default=800,help='较强暂定候选的MF严格下限，默认>800')
    app.add_argument('--rmf-min',type=score_threshold,default=800,help='纳入候选的RMF严格下限，默认>800')
    app.add_argument('--column',default='',help='采集色谱柱，如HP-5ms；未指定时不宣称外部RI条件已验证')
    app.add_argument('--blank-sample',action='append',default=[],help='峰表中CK的完整sample名称，可重复')
    app.add_argument('--evidence-json',type=absolute,help='增加已核查的外部RI/文献JSON；见文档示例')
    app.add_argument('--no-bundled-evidence',action='store_true',help='不加载随工具附带的公开文献记录')
    app.add_argument('--keep-peak',action='append',default=[],help='另保留此峰RMF过阈值的命中，仍标注RI或MF问题')


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
        print(f"{label} {original.name}：{len(scans)}次扫描，{len(p)}个峰。", flush=True)
    write_peak_csv(args.out / "peaks.csv", peaks)
    write_msp(args.out / "queries.msp", queries)
    return peaks, queries, native_scans, manifest


def existing_hits(args, msp, queries):
    if args.nist_dir:
        print(f"向 NIST MS Search 提交{len(queries)}张谱图。", flush=True)
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
    reference_index=[]
    if args.adams:
        print("读取Adams完整条目索引，并提取RI窗口内参考谱。", flush=True)
        references, ref_count = pdf_entries(source, ri, args.ri_window, reference_index)
    elif args.reference_json:
        references = local_reference_json(source)
        ref_count = len(references)
        reference_index=[{k:v for k,v in r.items() if k!='ions'} for r in references]
    else:
        references=[];ref_count=0
    matches, fragments, diagnostic = window_matches(
        peaks, references, queries, links, hits, anchors, args.ri_window, native_scans,
        (args.mass_min, args.mass_max),
    )
    pairs = replicate_pairs([p for p in peaks if p['sample'] not in args.blank_sample], queries, links)
    workbook = args.out / "report.xlsx"
    metadata = {
        "gcms-ri version": pkg_version(), "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "python": platform.python_version(), "ri_window": args.ri_window,
        "reference": str(source) if source else None, "reference_sha256": sha256(source) if source else None, "reference_index_count": ref_count,
        "reference_kind": "Adams" if args.adams else "local",
        "column": args.column, "blank_samples": args.blank_sample,
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
    for name,records in [('reference-index.json',reference_index),('reference-spectra.json',references)]:
        (args.out/name).write_text(json.dumps(records,ensure_ascii=False,indent=2),encoding='utf-8')
        metadata[name+'_sha256']=sha256(args.out/name)
    emit_summary(args,peaks,queries,links,hits,anchors,reference_index,references,metadata)
    (args.out / "run.json").write_text(json.dumps(metadata, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"{len(peaks)}个峰，{len(matches)}条RI窗口记录，"
          f"{sum(m['cosine'] is not None for m in matches)}条已比较参考谱。")
    print(workbook)


def emit_summary(args,peaks,queries,links,hits,anchors,index,references,metadata):
    evidence=load_evidence(args.evidence_json,not args.no_bundled_evidence)
    payload=summarize(peaks,queries,links,hits,anchors,index,references,evidence,
        source_kind=metadata.get('reference_kind','local'),column=args.column,blank_samples=args.blank_sample,
        mf_min=args.mf_min,rmf_min=args.rmf_min,window=args.ri_window,
        mass_range=(args.mass_min,args.mass_max),keep_peaks=args.keep_peak)
    metadata['candidate_settings']=payload['settings']
    metadata['candidate_counts']=payload['counts']
    metadata['nist_exported_hit_counts']=sorted({len(h) for h in hits.values()})
    if args.evidence_json:metadata['evidence_json_sha256']=sha256(args.evidence_json)
    evidence_path=args.out/'evidence-used.json'
    evidence_path.write_text(json.dumps(evidence,ensure_ascii=False,indent=2),encoding='utf-8')
    metadata['evidence_used_sha256']=sha256(evidence_path)
    write_reports(args.out,payload,metadata)
    print(f"中文汇总：{payload['counts']['candidate_rows']}条候选记录、{payload['counts']['distinct_peaks']}个峰。")
    print(args.out/'候选汇总.xlsx')


def summarize_existing(args):
    source=args.run_dir
    peaks=peak_table(source/'peaks.csv');queries=parse_msp(source/'queries.msp')
    links=link_spectra(peaks,queries);anchors=load_anchors(source/'alkanes.csv')
    report=source/'nist/results.txt'
    hits=parse_report(report) if report.is_file() else {}
    if hits and hits.keys()!=queries.keys():
        raise ValueError('旧NIST报告与MSP查询名不一致，不能汇总。')
    old=json.loads((source/'run.json').read_text(encoding='utf-8')) if (source/'run.json').is_file() else {}
    references=[]
    if args.reference_index:
        index=json.loads(args.reference_index.read_text(encoding='utf-8-sig'))
        kind=args.index_kind
    elif (source/'reference-index.json').is_file():
        index=json.loads((source/'reference-index.json').read_text(encoding='utf-8'))
        if (source/'reference-spectra.json').is_file():
            references=json.loads((source/'reference-spectra.json').read_text(encoding='utf-8'))
        kind=old.get('reference_kind','local')
    elif old.get('reference') and Path(old['reference']).is_file():
        path=Path(old['reference']);index=[]
        if path.suffix.lower()=='.pdf':
            print('旧版本未缓存完整索引，正在读取本地Adams PDF。',flush=True)
            references,_=pdf_entries(path,[interpolate(p['rt_min'],anchors) for p in peaks],args.ri_window,index)
            kind='Adams'
        else:
            references=local_reference_json(path);index=references;kind='local'
    else:
        index=[];kind='local'
    if not isinstance(index,list):raise ValueError('参考索引须为JSON数组。')
    for row in index:
        if not isinstance(row,dict) or not row.get('name') or not math.isfinite(float(row['ai'])) or float(row['ai'])<=0:
            raise ValueError('参考索引每项须有name和正数ai。')
    if not args.blank_sample:args.blank_sample=old.get('blank_samples',[])
    if not args.column:args.column=old.get('column','')
    if not args.evidence_json and not args.no_bundled_evidence and (source/'evidence-used.json').is_file():
        # Preserve manually reviewed evidence when rerunning a cached analysis.
        args.evidence_json=source/'evidence-used.json'
        args.no_bundled_evidence=True
    load_evidence(args.evidence_json,not args.no_bundled_evidence)
    create_output(args.out)
    for name in ('peaks.csv','queries.msp','alkanes.csv'):
        shutil.copyfile(source/name,args.out/name)
    if report.is_file():
        (args.out/'nist').mkdir();shutil.copyfile(report,args.out/'nist/results.txt')
    for name,records in [('reference-index.json',index),('reference-spectra.json',references)]:
        (args.out/name).write_text(json.dumps(records,ensure_ascii=False,indent=2),encoding='utf-8')
    metadata={'gcms-ri version':pkg_version(),'timestamp_utc':datetime.now(timezone.utc).isoformat(),
        'reference_kind':kind,'reference_index_count':len(index),'source_run':str(source),
        'source_run_sha256':sha256(source/'run.json') if (source/'run.json').is_file() else None,
        'query_sha256':sha256(source/'queries.msp'),'calibration_sha256':sha256(source/'alkanes.csv'),
        'nist_report_sha256':sha256(report) if report.is_file() else None,
        'reference-index.json_sha256':sha256(args.out/'reference-index.json'),
        'reference-spectra.json_sha256':sha256(args.out/'reference-spectra.json'),
        'column':args.column,'blank_samples':args.blank_sample,
        'operation':'复用已有NIST报告，只重做候选汇总；未重新检索NIST'}
    emit_summary(args,peaks,queries,links,hits,anchors,index,references,metadata)
    (args.out/'run.json').write_text(json.dumps(metadata,ensure_ascii=False,indent=2),encoding='utf-8')


def execute(args):
    if args.command == 'literature':
        from .webbook import lookup
        create_output(args.out)
        results=lookup(args.cas,args.out,args.timeout)
        if any(r['status'].startswith('检索失败') for r in results):
            raise ValueError('部分公开RI查询失败，已保留日志和成功记录。失败不表示无文献。')
        return
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
    if args.command == 'summarize':
        summarize_existing(args)
        return
    # Validate evidence before running costly extraction/NIST stages.
    load_evidence(args.evidence_json,not args.no_bundled_evidence)
    anchors = load_anchors(args.alkanes)
    if args.command == "run":
        if args.skip_nist and (args.nist_dir or args.nist_report):
            raise ValueError("Choose NIST input or --skip-nist.")
        if not (args.skip_nist or args.nist_dir or args.nist_report):
            raise ValueError("Use --nist-dir, --nist-report, or --skip-nist explicitly for native runs.")
        if args.sample_dir and not args.sample_dir.is_dir():
            raise ValueError(f"Sample directory not found: {args.sample_dir}")
        paths = args.sample or sorted(p for p in args.sample_dir.iterdir() if p.is_dir() and p.suffix.lower() == ".d")
        paths=list(dict.fromkeys([*paths,*args.blank]))
        args.blank_sample=list(dict.fromkeys([*args.blank_sample,*(p.stem for p in args.blank)]))
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
    # Redirected Windows consoles can use a code page that cannot print some
    # sample paths. Keep the report on disk even when that path needs escaping.
    for stream in (sys.stdout, sys.stderr):
        if isinstance(stream, TextIOWrapper):
            stream.reconfigure(encoding=stream.encoding if stream.isatty() else 'utf-8',errors="backslashreplace")
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
