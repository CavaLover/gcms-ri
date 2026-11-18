"""Portable Chinese candidate exports; no proprietary spreadsheet runtime."""

import csv
import html
import json

import xlsxwriter


RULES = [
 'MF、RMF为谱图匹配分数，不是鉴定准确率；所有命名均为暂定候选。',
 '默认严格使用RMF>800与|实验RI−初筛参考RI|<20；MF≤800或只有RAW支持另标待复核。原详细RI候选表保留含边界的窗口，用途与本汇总不同。',
 'Adams精油库未检出该CAS不是直接排除理由。只在本次完整解析索引范围内描述覆盖情况，不代表所有精油库均无该化合物。',
 '实验RI由相邻、连续碳数的正构烷烃时间线性插值，标定范围外不外推。RI差值=实验值−参考值。',
 '高MF/RMF但实验RI或参考RI缺失的记录仍保留为待复核；不能把正构烷烃定义RI当成本峰的实验RI。',
 'NIST报告的库RI未注明固定相及实测/预测属性，单凭数值接近不足以完成外部RI验证。核查固定相类别、升温模式和文献出处。',
 '相近固定相的程序升温参考RI全部展示；同时存在支持与冲突时标记冲突，不挑选最接近实测值的外部记录。',
 '主要离子为采集范围内最强8项，按整张输入谱的基峰归一至100%；这些信号不自动等于结构诊断离子。',
 'BG为局部扣背景谱，RAW为原谱；均不是CK扣除或共洗脱解卷积。导入MSP的处理方法由导出来源决定。',
 'CK必须显式指定样品名；匹配使用RT差≤0.10 min、同处理模式谱余弦≥0.80。未匹配不等于空白不存在。峰顶信号比不是浓度比。CK不参与生物学重复配对。',
 '保留完整检索名中的(E)/(Z)等构型描述，但不会因此宣称已排除异构体。同峰多个CAS为竞争候选，候选记录数不等于成分数。',
 '本工具只处理实际导出的NIST命中；若每谱只导出前两项，则不能声称穷尽异构体。文献出现不证明本批样品身份或生物学功能。',
 '文献记录注明已读原文、经数据库追溯或仅数据库条目。未提供或未检索到资料时明确标注，不生成不存在的引用。',
]

HEADERS = ['峰号','样品','保留时间RT（min）','完整暂定候选名称','主要离子m/z（相对强度%）',
           'MF','RMF','实验RI','初筛参考RI','RI差值','初筛RI来源','精油库检索备注','复核状态','CK核查',
           '外部RI与文献编号','外部证据与限制','CAS','采用谱/排名','其他竞争候选','备注']


def ions_text(ions):
    return '；'.join(f"{i['mz']} ({i['percent']:.1f})" if i['percent'] is not None else f"{i['mz']} (未提取或范围外)" for i in ions)


def hit_text(hit):
    return f"{hit['name']} [{hit['cas']}] MF/RMF {hit['mf']}/{hit['rmf']}；库RI {hit['library_ri']}"


def candidate_values(row):
    blank = row['blank']['status']
    if row['blank']['best']:
        best=row['blank']['best']
        blank += f"：{best['blank_peak']}；{best['mode']}余弦{best['cosine']:.3f}；峰顶信号比{best['signal_ratio']:.2f}"
    refs, notes = [], []
    for e in row['external_evidence']:
        delta=f"{e['delta_ri']:+.1f}" if e['delta_ri'] is not None else '未计算'
        refs.append(f"{e['id']}：RI {e['ri'] if e['ri'] is not None else '未提供'}；差值{delta}；{e['column']}")
        notes.append(f"{e['id']}：{e['assessment']}；{e['verification']}；{e['note']}")
    return [row['peak_id'],row['sample'],row['rt_min'],row['name'],ions_text(row['sample_ions']),
            row['mf'],row['rmf'],row['ri'],row['reference_ri'],row['delta_ri'],row['reference_source'],
            row['coverage'],row['grade'],blank,'\n'.join(refs) or '未核实对应外部RI',
            '\n'.join(notes) or row['literature_status'],row['cas'],f"{row['mode']} 第{row['rank']}项",
            '\n'.join(f"{h['mode']} #{h['rank']}：{hit_text(h)}" for h in row['competing_hits']),
            '\n'.join(row['caveats'])]


def write_reports(out, payload, metadata):
    payload['metadata'] = metadata
    payload['rules'] = RULES
    (out/'candidates.json').write_text(json.dumps(payload,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    rows=[candidate_values(r) for r in payload['candidates']]
    with (out/'候选汇总.csv').open('w',encoding='utf-8-sig',newline='') as f:
        writer=csv.writer(f);writer.writerow(HEADERS)
        # Prevent formula injection when a CSV is opened in a spreadsheet.
        writer.writerows(["'"+v if isinstance(v,str) and v.startswith(('=','+','-','@')) else v for v in row] for row in rows)
    used={e['id']:e for r in payload['candidates'] for e in r['external_evidence']}
    sources=[[e['id'],e['cas'],e.get('name',''),e['ri'],e['column'],
              {'programmed':'程序升温','isothermal':'恒温','unknown':'未核实'}[e['temperature_mode']],
              e['citation'],e['verification'],e['note'],e['url'],e['source_url']] for e in used.values()]
    source_headers=['编号','CAS','化合物','参考RI','色谱柱','升温模式','论文或资料','核查范围','说明','论文/资料链接','数据来源链接']
    checks=[]
    for r in payload['candidates']:
        diagnostic='\n'.join(e['id']+'：'+ions_text(e['measured_ions']) for e in r['external_evidence'] if e['measured_ions'])
        raw='\n'.join(hit_text(h) for h in r['raw_bg_hits']['RAW'])
        bg='\n'.join(hit_text(h) for h in r['raw_bg_hits']['BG'])
        checks.append([r['peak_id'],r['name'],ions_text(r['sample_ions']),
                       ions_text(r['reference_ions']) or '未导出对应参考谱强度',diagnostic or '未核实诊断离子',raw,bg])
    partial=out/'候选汇总.partial.xlsx'
    with xlsxwriter.Workbook(partial,{'strings_to_formulas':False,'strings_to_urls':False}) as wb:
        title=wb.add_format({'bold':True,'font_color':'white','bg_color':'#254D6B','text_wrap':True,'valign':'vcenter'})
        text=wb.add_format({'text_wrap':True,'valign':'vcenter','font_name':'Microsoft YaHei','font_size':10})
        decimal=wb.add_format({'num_format':'0.0','valign':'vcenter'})
        rt=wb.add_format({'num_format':'0.0000','valign':'vcenter'})
        amber=wb.add_format({'bg_color':'#FFF3DB','text_wrap':True,'valign':'vcenter'})
        def table(name,headers,records,widths):
            s=wb.add_worksheet(name);s.hide_gridlines(2);s.freeze_panes(1,2)
            s.set_row(0,34);s.write_row(0,0,headers,title)
            for j,width in enumerate(widths):s.set_column(j,j,width,text)
            for i,row in enumerate(records,1):
                s.set_row(i,110);s.write_row(i,0,row,text)
            s.autofilter(0,0,len(records),len(headers)-1)
            return s
        sheet=table('候选总表',HEADERS,rows,[16,24,16,48,38,9,9,12,12,12,30,44,28,40,56,85,16,16,74,70])
        for i,row in enumerate(rows,1):
            sheet.write_number(i,2,row[2],rt)
            for j in (7,8,9):
                if row[j] is not None:sheet.write_number(i,j,row[j],decimal)
            if '待复核' in row[12]:sheet.write_string(i,12,row[12],amber)
        table('碎片与命中', ['峰号','完整候选名','采用谱离子m/z（%）','参考谱对应离子m/z（%）','文献所列离子在样品中的强度（%）','RAW命中','BG命中'],checks,[16,48,45,45,65,75,75])
        table('文献索引',source_headers,sources,[23,16,38,12,18,16,80,42,85,65,75])
        info=[['候选记录数',payload['counts']['candidate_rows']],['覆盖峰数',payload['counts']['distinct_peaks']],
              ['版本',metadata.get('gcms-ri version','source')],['采集柱',payload['settings']['column'] or '未指定'],
              ['CK样品','；'.join(payload['settings']['blank_samples']) or '未指定'],
              ['阈值',f"MF>{payload['settings']['mf_min']}；RMF>{payload['settings']['rmf_min']}；|ΔRI|<{payload['settings']['ri_window']}"],
              *[[f'规则{i}',rule] for i,rule in enumerate(RULES,1)]]
        table('筛选说明',['项目','说明'],info,[28,145])
    partial.replace(out/'候选汇总.xlsx')
    esc=lambda v:html.escape(str(v if v is not None else '—'))
    body='''<!doctype html><html lang="zh"><meta charset="utf-8"><title>GC–MS候选复核</title>
<style>body{font:14px 'Microsoft YaHei',sans-serif;margin:24px;color:#17364b}table{border-collapse:collapse;width:100%;font-size:13px}th,td{padding:9px;border:1px solid #d7e0e6;white-space:pre-line;vertical-align:top}th{background:#254d6b;color:white;position:sticky;top:0}tr:nth-child(even){background:#f1f6f8}input{padding:8px;width:420px;margin:14px 0}.scroll{overflow:auto}td{min-width:90px}td:nth-child(4),td:nth-child(5),td:nth-child(16){min-width:280px}p{line-height:1.7;max-width:1150px}</style><h1>GC–MS候选与文献复核</h1>'''
    body+=f"<p>{len(rows)}条暂定候选记录、{payload['counts']['distinct_peaks']}个峰。候选数不等于已鉴定成分数；MF/RMF不是准确率。</p>"
    body+='<input id="q" aria-label="筛选候选" placeholder="输入样品、化合物、CAS或状态"><div class="scroll"><table id="candidates"><thead><tr>'
    body+=''.join('<th>'+esc(h)+'</th>' for h in HEADERS)+'</tr></thead><tbody>'
    for row in rows:
        body+='<tr>'+''.join('<td>'+esc(round(v,4) if isinstance(v,float) else v)+'</td>' for v in row)+'</tr>'
    body+='</tbody></table></div><h2>文献与资料</h2>'
    for e in used.values():
        body+='<p><b>'+esc(e['id'])+'</b> '+esc(e['citation'])+'<br>'+esc(e['verification'])+'；'+esc(e['note'])+f'<br><a href="{esc(e["url"])}">论文或资料</a> · <a href="{esc(e["source_url"])}">数据出处</a></p>'
    body+='<h2>筛选说明</h2>'+''.join('<p>'+esc(x)+'</p>' for x in RULES)
    body+='''<script>document.getElementById('q').addEventListener('input',e=>{let q=e.target.value.toLowerCase();document.querySelectorAll('#candidates tbody tr').forEach(r=>r.hidden=!r.textContent.toLowerCase().includes(q));});</script></html>'''
    (out/'候选汇总.html').write_text(body,encoding='utf-8')
