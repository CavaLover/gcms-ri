"""Candidate inclusion independent of essential-oil-library coverage."""

from collections import Counter, defaultdict
import math

from .evidence import assess
from .review import CAS
from .spectra import cosine, interpolate, relative


def number(value):
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None
    return result if math.isfinite(result) and result > 0 else None


def top_ions(values, mass_range, count=8):
    values = relative(values)
    lower, upper = mass_range
    masses = sorted((m for m in range(lower, upper + 1) if values[m] > 0), key=lambda m: (-values[m], m))[:count]
    return [{'mz': m, 'percent': float(values[m])} for m in masses]


def blank_matches(peaks, queries, links, blank_samples):
    available = {p['sample'] for p in peaks}
    if set(blank_samples) - available:
        raise ValueError(f'空白样品名未出现在峰表中：{sorted(set(blank_samples) - available)}')
    blanks = [p for p in peaks if p['sample'] in blank_samples]
    result = {}
    for peak in peaks:
        if peak['sample'] in blank_samples:
            continue
        modes = links[peak['peak_id']]
        candidates = []
        for blank in blanks:
            if abs(peak['rt_min'] - blank['rt_min']) > .10:
                continue
            blank_modes = links[blank['peak_id']]
            # Compare the same processing mode; never mix RAW with BG.
            mode = 'BG' if 'BG' in modes and 'BG' in blank_modes else 'RAW'
            if mode not in modes or mode not in blank_modes:
                continue
            a, b = queries[modes[mode]], queries[blank_modes[mode]]
            score = cosine(a, b)
            if score is not None and score >= .80:
                candidates.append(dict(blank_peak=blank['peak_id'], cosine=score, mode=mode,
                                       signal_ratio=float(a.sum()/b.sum()) if b.sum() else None))
        best = max(candidates, key=lambda c:c['cosine']) if candidates else None
        result[peak['peak_id']] = {'peak_id': peak['peak_id'], 'matches': candidates,
            'status': ('未指定CK' if not blank_samples else '匹配到CK检出峰' if best else
                       '未匹配到CK检出峰（不等于空白不存在）'), 'best': best}
    return result


def summarize(peaks, queries, links, hits, anchors, reference_index, references, evidence,
              *, source_kind='Adams', column='', blank_samples=(), mf_min=800, rmf_min=800,
              window=20, mass_range=(38,600), keep_peaks=()):
    by_cas = defaultdict(list)
    for entry in reference_index:
        if CAS.fullmatch(entry.get('cas', '')):
            by_cas[entry['cas']].append(entry)
    if set(keep_peaks) - {p['peak_id'] for p in peaks}:
        raise ValueError('指定保留的峰号不在峰表中。')
    blanks = blank_matches(peaks, queries, links, blank_samples)
    output, audit = [], []
    for peak in peaks:
        if peak['sample'] in blank_samples:
            continue
        pid = peak['peak_id']
        measured = interpolate(peak['rt_min'], anchors)
        modes = links[pid]
        grouped = defaultdict(list)
        for mode in ('BG','RAW'):
            for hit in hits.get(modes.get(mode), []):
                grouped[hit['cas']].append((mode, hit))
        for cas, choices in grouped.items():
            if not CAS.fullmatch(cas):
                audit.append(dict(peak_id=pid,cas=cas,included=False,reason='无可用CAS，未自动定名'))
                continue
            known = by_cas.get(cas, [])
            matching = [r for r in known if measured is not None and abs(measured-float(r['ai'])) < window]
            matching.sort(key=lambda r:abs(measured-float(r['ai'])))
            options = []
            for mode, hit in choices:
                ref = matching[0] if matching else None
                value = float(ref['ai']) if ref else number(hit.get('library_ri'))
                delta = measured-value if measured is not None and value is not None else None
                high_reverse = hit['rmf'] > rmf_min
                high_forward = hit['mf'] > mf_min
                passes_ri = delta is not None and abs(delta) < window
                unvalidated = high_forward and (measured is None or value is None)
                if high_reverse and (passes_ri or unvalidated or pid in keep_peaks):
                    options.append((mode,hit,ref,value,delta,passes_ri))
            if not options:
                audit.append(dict(peak_id=pid,cas=cas,included=False,reason='未通过RMF与RI初筛；且不属高MF但RI缺失或手动保留记录'))
                continue
            # Prefer background-corrected matches, then their exported rank.
            mode,hit,ref,value,delta,passes_ri = min(options,key=lambda x:(x[0]!='BG',x[1]['rank']))
            external = assess(evidence,cas,measured,column,window)
            compatible = [e for e in external if e['comparable'] and e['delta_ri'] is not None]
            disagreements = [e for e in compatible if abs(e['delta_ri']) >= window]
            supports = [e for e in compatible if abs(e['delta_ri']) < window]
            caveats = ['暂定候选；未用同条件标准品确证，未穷尽同分异构体。']
            if source_kind == 'Adams':
                coverage = ('精油库检出对应CAS' if known else
                            '精油库里面搜索不到：本次Adams索引未检出该CAS对应条目')
            elif reference_index:
                coverage = ('本地参考表检出对应CAS；未据此判断Adams覆盖情况' if known else
                            '本地参考表未检出对应CAS；Adams精油库未核查')
            else:
                coverage = '未提供Adams或参考索引，精油库覆盖情况未核查'
            grade = '较强暂定候选'
            if measured is None:
                grade = '待复核：实验RI无法内插'
            elif value is None:
                grade = '待复核：初筛参考RI缺失'
            elif not passes_ri:
                grade = '待复核：手动保留但RI不一致'
            elif disagreements:
                grade = '待复核：外部RI存在冲突'
            elif known and not matching:
                grade = '待复核：本地参考AI不一致'
            elif hit['mf'] <= mf_min or mode == 'RAW':
                grade = '待复核：MF不足或仅RAW支持'
            elif not ref and not supports:
                grade = '待复核：外部RI条件未验证'
            if not ref:
                caveats.append('NIST报告库RI未附固定相及实测/预测属性；数值接近只作初筛。')
            if known and not matching:
                caveats.append('参考库有该CAS，但AI未通过窗口；不能备注为精油库搜不到。')
            if disagreements and supports:
                caveats.append('外部可比RI同时有支持和冲突，全部列出，不选择最接近值掩盖冲突。')
            if measured is None:
                caveats.append('标定范围外不外推，正构烷烃定义RI不能代替本峰的实验RI。')
            related = []
            for other_mode in ('BG','RAW'):
                for h in hits.get(modes.get(other_mode),[]):
                    if h['cas'] != cas:
                        related.append(dict(mode=other_mode, **h))
            spectrum = queries[modes[mode]]
            sample_ions = top_ions(spectrum,mass_range)
            reference_ions = []
            if ref:
                spectra = [r for r in references if r.get('cas')==cas and float(r['ai'])==float(ref['ai']) and r.get('ions')]
                if spectra:
                    record = spectra[0]; ions={int(m):float(a) for m,a in record['ions'].items()}; maximum=max(ions.values())
                    reference_ions=[dict(mz=i['mz'],percent=100*ions[i['mz']]/maximum if i['mz'] in ions else None) for i in sample_ions]
                    caveats.append('Adams PDF提取的参考强度为近似值；未提取到不等于不存在。' if source_kind=='Adams' else '参考强度来自提供的本地参考谱。')
            v=relative(spectrum)
            for e in external:
                e['measured_ions']=[dict(mz=m,percent=float(v[m]) if mass_range[0]<=m<=mass_range[1] else None) for m in e['ions']]
            evidence_note='；'.join(dict.fromkeys(e['note'] for e in external if e['note']))
            raw_hits = {m:[h for mode_,h in choices if mode_==m] for m in ('RAW','BG')}
            row=dict(peak_id=pid,sample=peak['sample'],rt_min=peak['rt_min'],ri=measured,
                     cas=cas,name=hit['name'],mode=mode,query=modes[mode],mf=hit['mf'],rmf=hit['rmf'],rank=hit['rank'],
                     library=hit['library'],library_id=hit['library_id'],reference_ri=value,delta_ri=delta,
                     reference_source=('Adams AI' if source_kind=='Adams' else '本地参考AI') if ref else 'NIST报告库RI（条件未注明）',
                     reference_page=ref.get('pdf_page') if ref else None,coverage=coverage,grade=grade,
                     sample_ions=sample_ions,reference_ions=reference_ions,external_evidence=external,
                     literature_status='已列对应资料及核查范围' if external else '本轮未检出已核实的对应文献记录',
                     blank=blanks[pid],competing_hits=related,raw_bg_hits=raw_hits,caveats=caveats,
                     evidence_note=evidence_note,manual_keep=pid in keep_peaks,passes_initial_ri=passes_ri)
            output.append(row)
            audit.append(dict(peak_id=pid,cas=cas,included=True,reason=grade))
    # Same peak, different CAS: retain each as an alternative, never count them as separate detections.
    counts=Counter(r['peak_id'] for r in output)
    for row in output:
        if counts[row['peak_id']] > 1:
            row['caveats'].append(f"同峰保留{counts[row['peak_id']]}个CAS候选，不能计为多个已检出成分。")
            if row['grade']=='较强暂定候选' and row['rank']>1:
                row['grade']='同峰竞争候选'
    output.sort(key=lambda r:(r['sample'],r['rt_min'],r['rank'],r['cas']))
    return dict(schema_version=1,candidates=output,audit=audit,blank_matches=list(blanks.values()),
                settings=dict(mf_min=mf_min,rmf_min=rmf_min,ri_window=window,mass_range=list(mass_range),
                              column=column,blank_samples=list(blank_samples),keep_peaks=list(keep_peaks),
                              source_kind=source_kind,strict_thresholds=True),
                counts=dict(candidate_rows=len(output),distinct_peaks=len(counts),
                            grades=dict(Counter(r['grade'] for r in output))))
