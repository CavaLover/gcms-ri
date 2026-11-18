"""Optional public NIST RI lookup. No spectra or sample data are uploaded."""

from datetime import datetime, timezone
from html.parser import HTMLParser
import json
import re
import time
from urllib.error import URLError
from urllib.request import Request, urlopen

from .review import CAS, sha256


class RetentionTables(HTMLParser):
    def __init__(self):
        super().__init__(convert_charrefs=True)
        self.tables=[]; self.heading=''; self.in_heading=False; self.heading_text=[]
        self.table=None; self.row=None; self.cell=None

    def handle_starttag(self,tag,attrs):
        if tag in ('h2','h3','h4'):
            self.in_heading=True;self.heading_text=[]
        if tag=='table':
            self.table={'section':self.heading,'rows':[]}
        if self.table is not None and tag=='tr':self.row=[]
        if self.row is not None and tag in ('td','th'):self.cell=[]

    def handle_data(self,data):
        if self.in_heading:self.heading_text.append(data)
        if self.cell is not None:self.cell.append(data)

    def handle_endtag(self,tag):
        if tag in ('h2','h3','h4') and self.in_heading:
            self.heading=' '.join(''.join(self.heading_text).split());self.in_heading=False
        if tag in ('td','th') and self.cell is not None:
            self.row.append(' '.join(''.join(self.cell).split()));self.cell=None
        if tag=='tr' and self.table is not None and self.row is not None:
            self.table['rows'].append(self.row);self.row=None
        if tag=='table' and self.table is not None:
            self.tables.append(self.table);self.table=None


def parse_webbook(text,cas,url):
    parser=RetentionTables();parser.feed(text)
    # Avoid accepting a different CAS returned after a failed identifier lookup.
    if not re.search(r'CAS Registry Number\s*:\s*(?:</?[^>]+>\s*)*'+re.escape(cas)+r'\b',text,re.I):
        raise ValueError('返回页面未验证到请求的CAS，未导入其RI。')
    records=[]
    for table in parser.tables:
        rows=table['rows']
        if not rows:continue
        headers=rows[0]
        if 'Active phase' not in headers or 'I' not in headers:continue
        section=table['section']
        mode='isothermal' if 'isothermal' in section.lower() else 'programmed' if ('temperature ramp' in section.lower() or 'temperature program' in section.lower()) else 'unknown'
        for row in rows[1:]:
            if len(row)!=len(headers):continue
            fields=dict(zip(headers,row))
            value=fields['I']
            if not re.fullmatch(r'\d+(?:\.\d*)?',value):continue
            reference=fields.get('Reference','NIST Chemistry WebBook')
            records.append(dict(id=f'webbook-{cas}-{len(records)+1}',cas=cas,ri=float(value),
                column=fields['Active phase'],temperature_mode=mode,citation=reference,url=url,source_url=url,
                verification='已提取NIST公开RI表，未核读原始论文全文',
                note=f"{section}；{fields.get('Comment','')}；数据库记录不等于样品确证。",ions=[]))
    return records


def lookup(cases,out,timeout=20):
    cases=list(dict.fromkeys(cases))
    if not cases or any(not CAS.fullmatch(cas) for cas in cases):
        raise ValueError('请提供有效CAS。')
    evidence=[]; log=[]
    for i,cas in enumerate(cases):
        if i:time.sleep(1)
        url='https://webbook.nist.gov/cgi/cbook.cgi?ID=C'+cas.replace('-','')+'&Mask=2000'
        entry=dict(cas=cas,url=url,checked=datetime.now(timezone.utc).isoformat())
        print(f'查询公开NIST RI资料：{cas}',flush=True)
        try:
            request=Request(url,headers={'User-Agent':'gcms-ri/0.2 (public retention-index research)'})
            with urlopen(request,timeout=timeout) as response:
                raw=response.read(8_000_001)
            if len(raw)>8_000_000:raise ValueError('网页超出8 MB，未处理。')
            path=out/f'{cas}.html';path.write_bytes(raw)
            records=parse_webbook(raw.decode('utf-8'),cas,url)
            evidence.extend(records)
            entry.update(status='已提取' if records else '页面已读取，但未提取到RI记录',records=len(records),sha256=sha256(path))
        except (URLError,OSError,ValueError,UnicodeError) as exc:
            entry.update(status='检索失败，不能解释为没有文献',error=str(exc),records=0)
        log.append(entry)
    (out/'evidence.json').write_text(json.dumps(evidence,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    (out/'lookup.json').write_text(json.dumps(log,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(f'提取{len(evidence)}条公开RI记录；详情见 {out / "lookup.json"}')
    return log
