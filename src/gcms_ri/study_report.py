"""Export study worksheets and the blank evidence used to screen each peak."""

from collections import defaultdict

import xlsxwriter

from .report import _table


HEADERS = ("TR","计算得RI","adams RI","adams检索结果","nist检索结果",
           "Ions","相对分子质量","Formule","样品","备注")
WIDTHS = (14,14,13,28,44,31,15,16,38,88)


def format_ions(ions):
    return " / ".join(f"{ion['mz']}({ion['percent']:.0f}%)" for ion in ions[:6])


def result_row(peak):
    nist = peak["nist"].get("BG") or peak["nist"].get("RAW") or []
    entry = peak["adams"]
    note = [peak["decision"]]
    if peak["candidate"]:
        note.append(f"综合暂定：{peak['candidate']}")
    if nist:
        first = nist[0]
        note.append(f"NIST 首选 MF/RMF {first['mf']}/{first['rmf']}")
    alkane = peak["alkane"]
    if alkane:
        disposition = "同峰支持" if alkane["supported"] else "分子离子证据不足"
        note.append(f"n-C{alkane['carbon']}，标样 TR {alkane['standard_rt_min']:.4f} min，"
                    f"分子离子 m/z {alkane['parent_mz']}，{disposition}")
    if entry:
        note.append(f"Adams 同 CAS 条目，PDF 页 {entry['pdf_page']}" if entry.get("pdf_page") is not None
                    else "Adams 同 CAS 条目")
    elif peak["adams_alternatives"]:
        values = " / ".join(str(ref["ai"]) for ref in peak["adams_alternatives"])
        note.append(f"Adams 同 CAS 有多个 RI：{values}")
    elif peak["cas"] and peak["adams_index_available"]:
        note.append("本次 Adams 索引未返回对应条目")
    elif peak["cas"]:
        note.append("未提供 Adams 索引")
    if peak["ri"] is None:
        note.append("超出校准边界")
    if len(peak["alternatives"]) == 1 and not (alkane and alkane["supported"]):
        note.append(peak["alternatives"][0]["grade"])
    if len(peak["alternatives"]) > 1:
        note.append("同峰竞争候选："+" / ".join(r["name"] for r in peak["alternatives"]))
    return (peak["rt_min"],peak["ri"],entry["ai"] if entry else None,entry["name"] if entry else None,
            nist[0]["name"] if nist else None,format_ions(peak["ions"]),peak["molecular_weight"],
            peak["formula"],f"{peak['sample']} [{peak['id']}]", "。".join(note)+"。")


def write_study(path, data):
    groups = defaultdict(list)
    members = defaultdict(list)
    for sample in data["samples"]:
        members[sample["group"]].append(sample)
        groups.setdefault(sample["group"],[])  # Include groups whose peaks were all screened out.
    for peak in data["peaks"]:
        if not peak["excluded"]:
            groups[peak["group"]].append(peak)
    calibration = data["calibration"]
    with xlsxwriter.Workbook(path, {"strings_to_formulas":False, "strings_to_urls":False}) as book:
        head = book.add_format(dict(bold=True,font_color="#FFFFFF",bg_color="#254D6B"))
        base = book.add_format(dict(text_wrap=True,valign="vcenter"))
        numeric = book.add_format(dict(num_format="0.0000",valign="vcenter"))
        blue = book.add_format(dict(bg_color="#BDDDF8",text_wrap=True,valign="vcenter"))
        blue_number = book.add_format(dict(bg_color="#BDDDF8",num_format="0.0000",valign="vcenter"))
        standard_rows = []
        for anchor in calibration["anchors"]:
            n, entry = anchor["carbon"], anchor["adams"]
            row = (anchor["rt_min"],anchor["ri"],entry["ai"] if entry else None,
                   entry["name"] if entry else None,None,None,14*n+2,f"C{n}H{2*n+2}",
                   f"n-C{n} 标样 · {calibration['file']}", "使用者实测校准点。RI 按碳数定义。")
            standard_rows.append((anchor["rt_min"],True,row))
        mapping = []
        for i,(group,peaks) in enumerate(groups.items(),1):
            sheet = book.add_worksheet(f"Sheet{i}")
            mapping.append((sheet.name,group,len(members[group]),len(peaks)))
            sheet.freeze_panes(1,0)
            sheet.hide_gridlines(2)
            for j,width in enumerate(WIDTHS):
                sheet.set_column(j,j,width)
            sheet.write_row(0,0,HEADERS,head)
            entries = standard_rows + [(p["rt_min"],False,result_row(p)) for p in peaks]
            for j,(_,is_standard,row) in enumerate(sorted(entries,key=lambda item:(item[0],not item[1])),1):
                if j >= 1_048_576:
                    raise ValueError("单工作表超出 Excel 行限制。请分组或缩小数据集。")
                for col,value in enumerate(row):
                    style = (blue_number if col == 0 else blue) if is_standard else (numeric if col == 0 else base)
                    sheet.write(j,col,value,style)
            sheet.autofilter(0,0,len(entries),len(HEADERS)-1)
        _table(book,"分组索引",("工作表","分组","样品数","保留峰数"),mapping,(15,31,15,15),head,{})
        homologs = []
        for group, peaks in groups.items():
            detections = defaultdict(set)
            weak = defaultdict(set)
            for peak in peaks:
                alkane = peak["alkane"]
                if alkane:
                    target = detections if alkane["supported"] else weak
                    target[alkane["carbon"]].add(peak["sample_id"])
            for n in sorted(detections.keys() | weak.keys()):
                homologs.append((group,n,len(members[group]),len(detections[n]),
                                 ", ".join(sorted(detections[n])),", ".join(sorted(weak[n]-detections[n]))))
        _table(book,"分组烷烃对照",("分组","碳数","分组样品数","有分子离子支持的样品数","对应样品","其他弱线索"),
               homologs,(26,14,18,29,33,33),head,{})
        _table(book,"样品与CK",("样品编号","样品","分组","运行目录","原始扫描可用","CK来源","CK范围","CK原始扫描可用"),
               ((s["id"],s["sample"],s["group"],s["run_dir"],s["has_scans"],b["source"],b["scope"],b["has_scans"])
                for s in data["samples"] for b in s["blanks"]),
               (13,28,24,58,18,65,17,22),head,{})
        _table(book,"样品间峰配对",("样品A","样品B","峰A","峰B","校正ΔTR","EI余弦","谱处理"),
               (tuple(p[k] for k in ("sample_a","sample_b","peak_a","peak_b","delta_rt_min","cosine","mode"))
                for p in data["pairs"]),(15,15,22,22,17,17,14),head,{})
        _table(book,"CK审计",("峰号","分组","TR","状态","离子","样品EIC面积","CK",
                             "CK峰号","CK TR","EI余弦","峰形r","CK EIC面积","面积比"),
               ((p["id"],p["group"],p["rt_min"],p["decision"],"/".join(map(str,p["diagnostic_ions"])),
                 p["sample_area"],*(b.get(k) for k in ("blank","peak_id","rt_min","cosine","shape_r","area","sample_blank_ratio")))
                for p in data["peaks"] for b in p["blank_matches"] or [{}]),
               (22,24,12,31,16,16,70,20,14,16,16,16,16),head,{})
        rules = [
            ("版本",data["version"]),
            ("标样来源",calibration["file"]),
            ("标样文件 SHA-256",calibration["sha256"]),
            ("蓝色行","标样实测 TR 和按碳数定义的 RI；与样品峰按 TR 穿插。"),
            ("保留指数","样品 RI 用输入的正构烷烃锚点线性内插，范围外留空。"),
            ("空白判断","保留时间和同处理 EI 谱先配对，有扫描时补查 EIC。CK 同峰或典型背景移出主表，证据见 CK 审计。"),
            ("EIC 面积","五扫描平滑，以峰顶 ±0.065 min 和局部线性基线计算。面积比需峰形相符，未扣除 CK、未换算浓度。"),
            ("历史 CK","代表性背景筛查，以使用者声明的适用范围解释。"),
            ("强度范围","谱表含原始与扣背景命中。试验样品和对应 CK 的扫描缺失时，仅做谱图复核。"),
            ("同峰定名","正构烷烃优先结合用户标样和分子离子同峰，结果为暂定。谱库分数及对应 Adams RI 独立记录。"),
            ("RT 偏移","设计中显式提供的偏移用于标样窗口和样品间配对。实验 RI 由未平移的 TR 计算。"),
        ]
        _table(book,"说明",("项目","内容"),rules,(25,118),head,{})
