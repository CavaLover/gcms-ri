"""Traceable public RI/literature records, with explicit method compatibility."""

from collections import defaultdict
from importlib.resources import files
import json
import math
import re

from .review import CAS


def load_evidence(path=None, bundled=True):
    records = []
    if bundled:
        records.extend(json.loads(files('gcms_ri').joinpath('data/literature.json').read_text(encoding='utf-8')))
    if path:
        supplied = json.loads(path.read_text(encoding='utf-8-sig'))
        if not isinstance(supplied, list):
            raise ValueError('文献 JSON 须为记录数组。')
        records.extend(supplied)
    seen = set()
    for record in records:
        if not isinstance(record, dict):
            raise ValueError('每条文献须为 JSON 对象。')
        for key in ('id', 'cas', 'citation', 'url', 'column', 'temperature_mode', 'verification'):
            if not isinstance(record.get(key), str) or not record[key].strip():
                raise ValueError(f'文献缺少文本字段：{key}')
        if not CAS.fullmatch(record['cas']):
            raise ValueError(f"文献 CAS 无效：{record['cas']}")
        if record['id'] in seen:
            raise ValueError(f"文献编号重复：{record['id']}；请使用独立编号。")
        seen.add(record['id'])
        if not record['url'].startswith(('https://', 'http://')):
            raise ValueError('文献 URL 须以 https:// 或 http:// 开头。')
        if record['temperature_mode'] not in ('programmed', 'isothermal', 'unknown'):
            raise ValueError('temperature_mode 须为 programmed、isothermal 或 unknown。')
        if record.get('ri') is not None:
            number = float(record['ri'])
            if not math.isfinite(number) or number <= 0:
                raise ValueError('参考 RI 须为有限正数或 null。')
            record['ri'] = number
        record.setdefault('ri', None)
        ions = record.get('ions', [])
        if not isinstance(ions, list) or any(type(m) is not int or not 1 <= m <= 1200 for m in ions):
            raise ValueError('文献 ions 须为 1–1200 的整数数组。')
        record.setdefault('ions', [])
        record.setdefault('note', '')
        record.setdefault('source_url', record['url'])
    return records


def phase_family(column):
    name = re.sub(r'[^a-z0-9]', '', column.lower()).replace('ultrainert', '').replace('ui', '')
    if name in ('hp5', 'hp5ms', 'db5', 'db5ms', 'vf5ms', 'rtx5', 'rtx5ms', 'bpx5', 'slb5ms', 'tr5ms', 'se54', 'ultra2'):
        return '5%-phenyl'
    if name in ('db1', 'db1ms', 'hp1', 'hp1ms', 'rtx1', 'ov1', 'ov101', 'se30', 'ultra1'):
        return 'methyl-silicone'
    if name in ('dbwax', 'hpinnowax', 'dbwaxetr', 'carbowax20m', 'supelcowax10'):
        return 'wax'
    return None


def assess(records, cas, measured_ri, column, window):
    selected = []
    family = phase_family(column)
    for record in records:
        if record['cas'] != cas:
            continue
        item = dict(record)
        item['delta_ri'] = measured_ri - record['ri'] if measured_ri is not None and record['ri'] is not None else None
        item['comparable'] = bool(family and family == phase_family(record['column']) and
                                  record['temperature_mode'] == 'programmed')
        if item['delta_ri'] is None:
            item['assessment'] = '未完成RI数值验证'
        elif not item['comparable']:
            item['assessment'] = '条件不可直接比较或尚未核实'
        elif abs(item['delta_ri']) < window:
            item['assessment'] = '相近固定相的程序升温RI支持'
        else:
            item['assessment'] = '相近固定相的程序升温RI不一致'
        selected.append(item)
    return selected
