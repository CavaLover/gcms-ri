"""Review declared sample groups and their blanks from cached GC–MS runs."""

from collections import Counter, defaultdict
from functools import lru_cache
from itertools import combinations
import json
from pathlib import Path

from .candidates import summarize, top_ions
from .chromatography import ScanProfile
from .nist import parse_report
from .review import CAS, sha256
from .screening import N_ALKANE_CAS, Blank, alkane_candidate, screen
from .spectra import SAFE_ID, cosine, interpolate, link_spectra, load_anchors, parse_msp, peak_table
from .study_report import write_study


def read_json(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


@lru_cache(maxsize=4)
def scan_profile(path):
    return ScanProfile(path) if path and path.is_file() else None


class Run:
    def __init__(self, path, anchors, window):
        self.path = path
        self.peaks = peak_table(path / "peaks.csv")
        self.queries = parse_msp(path / "queries.msp")
        links = link_spectra(self.peaks,self.queries)
        self.spectra = {pid:{k:self.queries[name] for k,name in modes.items()} for pid,modes in links.items()}
        self.report = path / "nist" / "results.txt"
        hits = parse_report(self.report) if self.report.is_file() else {}
        if hits and hits.keys() != self.queries.keys():
            raise ValueError(f"NIST 报告和 MSP 查询号不一致：{path}")
        self.hits = {pid:{k:hits.get(name,[]) for k,name in modes.items()} for pid,modes in links.items()}
        manifest = path / "run.json"
        self.metadata = read_json(manifest) if manifest.is_file() else {}
        self.calibration = path / "reference-index.json"
        records = read_json(self.calibration) if self.calibration.is_file() else []
        self.adams = records if self.metadata.get("reference_kind") == "Adams" else []
        self.by_sample = defaultdict(list)
        for peak in self.peaks:
            self.by_sample[peak["sample"]].append(peak)
        self.labels = {s["name"]:s["label"] for s in self.metadata.get("sample_sources", [])
                       if isinstance(s.get("name"),str) and SAFE_ID.fullmatch(s.get("label") or "")}
        self.candidates = defaultdict(list)
        review = summarize(self.peaks,self.queries,links,hits,anchors,records,[],[],
                           source_kind=self.metadata.get("reference_kind","local"),window=window)
        for row in review["candidates"]:
            self.candidates[row["peak_id"]].append(row)

    def get_peaks(self, sample):
        if sample not in self.by_sample and sample not in self.labels:
            raise ValueError(f"样品 {sample!r} 未在 {self.path} 的峰表或运行清单中找到。")
        return self.by_sample.get(sample,[])

    def get_profile(self, sample):
        label = self.labels.get(sample)
        return scan_profile(self.path / "scans" / f"{label}.scans.jsonl") if label else None


def adams_references(records):
    if not isinstance(records,list):
        raise ValueError("Adams 索引须为 JSON 数组。")
    by_cas = defaultdict(dict)
    for item in records:
        if not isinstance(item,dict) or not item.get("name") or not 0 < float(item["ai"]) < float("inf"):
            raise ValueError("Adams 索引条目须包含名称和有限正数 AI。")
        cas = item.get("cas") or ""
        if CAS.fullmatch(cas):
            by_cas[cas].setdefault(float(item["ai"]),item)
    return {cas:list(entries.values()) for cas,entries in by_cas.items()}


def unique_adams(index, cas):
    entries = index.get(cas,[])
    return entries[0] if len(entries) == 1 else None


def decide_candidate(run, peak, profile, anchors, rt_window, rt_offset, adams):
    pid = peak["peak_id"]
    modes = run.spectra[pid]
    spectrum = modes["BG"] if "BG" in modes else modes["RAW"]
    alkane = alkane_candidate(peak,spectrum,profile,anchors,rt_window,rt_offset)
    alternatives = run.candidates[pid]
    first = alternatives[0] if len(alternatives) == 1 else None
    cas = first["cas"] if first else None
    row = dict(candidate=None, cas=cas, molecular_weight=None, formula=None, alkane=alkane)
    if first:
        row["candidate"] = first["name"]
        evidence = next(h for h in first["raw_bg_hits"][first["mode"]] if h["rank"]==first["rank"])
        row["formula"] = evidence.get("formula")
        mass = evidence.get("mw") or ""
        if mass.isdigit():
            row["molecular_weight"] = int(mass)
    if alkane and alkane["supported"]:
        n = alkane["carbon"]
        row.update(candidate=f"正构烷烃 n-C{n}",cas=alkane["cas"],formula=f"C{n}H{2*n+2}",
                   molecular_weight=alkane["parent_mz"])
    row["adams"] = unique_adams(adams,row["cas"])
    row["adams_index_available"] = bool(adams)
    row["adams_alternatives"] = adams.get(row["cas"],[])
    row["alternatives"] = [{"cas":alt["cas"],"name":alt["name"],"grade":alt["grade"]} for alt in alternatives]
    return row


def pair_peaks(samples, spectra, window):
    pairs = []
    for (sa,left),(sb,right) in combinations(samples.items(),2):
        possible = []
        for a in left:
            for b in right:
                da = a["rt_min"]-a["rt_offset_min"]
                db = b["rt_min"]-b["rt_offset_min"]
                delta = db-da
                if abs(delta) > window:
                    continue
                x,y = spectra[a["id"]],spectra[b["id"]]
                mode = "BG" if "BG" in x and "BG" in y else "RAW"
                if mode not in x or mode not in y:
                    continue
                score = cosine(x[mode],y[mode])
                if score is not None and score >= .85:
                    possible.append((score,-abs(delta),a["id"],b["id"],delta,mode))
        used = set()
        for score,_,a_id,b_id,delta,mode in sorted(possible,reverse=True):
            if a_id in used or b_id in used:
                continue
            used.update((a_id,b_id))
            pairs.append(dict(sample_a=sa,sample_b=sb,peak_a=a_id,peak_b=b_id,
                              delta_rt_min=delta,cosine=score,mode=mode))
    return pairs


def study(args, create_output, version):
    anchors = load_anchors(args.alkanes)
    design = read_json(args.design)
    if not isinstance(design,dict) or not isinstance(design.get("samples"),list) or not design["samples"]:
        raise ValueError("设计文件须含非空 samples 数组。")
    runs = {}

    def run_for(name):
        path = (args.design.parent / Path(name).expanduser()).resolve()
        if path not in runs:
            runs[path] = Run(path,anchors,args.ri_window)
        return runs[path]

    selected = []
    for item in design["samples"]:
        if not isinstance(item,dict):
            raise ValueError("samples 每项须为对象。")
        run, name = run_for(item["run_dir"]), item["sample"]
        group = item["group"]
        if not isinstance(group,str) or not group.strip():
            raise ValueError("每个实验样品须有非空 group。")
        blanks = item["blanks"]
        if not isinstance(blanks,list) or not blanks:
            raise ValueError(f"{name} 须显式指定至少一份 CK。")
        checked = []
        for control in blanks:
            if not isinstance(control,dict):
                raise ValueError("CK 条目须为对象。")
            b_run = run_for(control["run_dir"]) if "run_dir" in control else run
            b_name = control["sample"]
            scope = control.get("scope","unspecified")
            if scope not in ("same_batch","historical","unspecified"):
                raise ValueError(f"未知 CK scope：{scope}")
            b_run.get_peaks(b_name)
            checked.append((b_run,b_name,scope))
        offset = float(item.get("rt_offset_min",0))
        if not -float("inf") < offset < float("inf"):
            raise ValueError("rt_offset_min 须为有限值。")
        if not run.get_peaks(name):
            raise ValueError(f"实验样品 {name!r} 没有可复核的峰。")
        selected.append((run,name,group,checked,offset))
    distinct = {(run.path,name) for run,name,*_ in selected}
    if len(distinct)!=len(selected):
        raise ValueError("同一运行目录中的实验样品不可重复。")
    if distinct & {(b_run.path,b_name) for _,_,_,checks,_ in selected for b_run,b_name,_ in checks}:
        raise ValueError("同一数据不能同时声明为实验样品与 CK。")
    records = read_json(args.adams_index) if args.adams_index else [
        ref for run in runs.values() for ref in run.adams
    ]
    adams = adams_references(records)
    experimental, spectra, sample_peaks, notes = [], {}, {}, []
    for i,(run,name,group,checks,offset) in enumerate(selected,1):
        sample_id = f"S{i}"
        profile = run.get_profile(name)
        blanks = [Blank(f"{b_run.path} :: {b_name}",scope,tuple(b_run.get_peaks(b_name)),
                        b_run.spectra,b_run.get_profile(b_name)) for b_run,b_name,scope in checks]
        kept = []
        for peak in run.get_peaks(name):
            pid = peak["peak_id"]
            modes = run.spectra[pid]
            hits = run.hits[pid].get("BG") or run.hits[pid].get("RAW") or []
            review = screen(peak,modes,hits,profile,blanks,args.rt_window)
            spectrum = modes["BG"] if "BG" in modes else modes["RAW"]
            row = dict(**peak,id=f"{sample_id}:{pid}",sample_id=sample_id,group=group,rt_offset_min=offset,
                       ri=interpolate(peak["rt_min"],anchors),ions=top_ions(spectrum,(1,1200)),
                       nist=run.hits[pid],**review)
            if not row["excluded"]:
                row.update(decide_candidate(run,peak,profile,anchors,args.rt_window,offset,adams))
                spectra[row["id"]] = modes
                kept.append(row)
            experimental.append(row)
        sample_peaks[sample_id] = kept
        notes.append(dict(id=sample_id,sample=name,group=group,run_dir=str(run.path),has_scans=profile is not None,
                          blanks=[dict(source=b.label,scope=b.scope,has_scans=b.profile is not None) for b in blanks]))
    payload = dict(schema_version=1,version=version,design_sha256=sha256(args.design),
                   calibration=dict(file=str(args.alkanes),sha256=sha256(args.alkanes),
                       anchors=[dict(carbon=n,rt_min=rt,ri=n*100,
                                     adams=unique_adams(adams,N_ALKANE_CAS.get(n))) for n,rt in anchors]),
                   settings=dict(rt_window=args.rt_window,ri_window=args.ri_window),
                   samples=notes, peaks=experimental, pairs=pair_peaks(sample_peaks,spectra,args.rt_window),
                   decisions=dict(Counter(row["decision"] for row in experimental)))
    create_output(args.out)
    write_study(args.out / "分组复核.xlsx",payload)
    (args.out / "study.json").write_text(json.dumps(payload,ensure_ascii=False,indent=2),encoding="utf-8")
    print(args.out / "分组复核.xlsx")
