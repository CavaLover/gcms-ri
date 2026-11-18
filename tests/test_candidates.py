import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
from urllib.error import HTTPError
from zipfile import ZipFile

import numpy as np

from gcms_ri.candidates import summarize
from gcms_ri.cli import main
from gcms_ri.evidence import load_evidence
from gcms_ri.nist import parse_report
from gcms_ri.spectra import link_spectra
from gcms_ri.webbook import lookup, parse_webbook


def hit(cas='6221-98-3',mf=870,rmf=890,ri='1451',rank=1,name='Synthetic candidate (E)-'):
    return dict(cas=cas,mf=mf,rmf=rmf,library_ri=ri,rank=rank,name=name,library='synthetic',library_id='1')


def reference(cas='6221-98-3',ri=1450,column='DB-5ms',mode='programmed',id='test'):
    return dict(id=id,cas=cas,ri=ri,column=column,temperature_mode=mode,citation='Synthetic citation',
                url='https://example.org/test',source_url='https://example.org/test',
                verification='合成测试，不是真实文献',note='',ions=[71,89,900])


class CandidateTest(unittest.TestCase):
    def setUp(self):
        self.peaks=[dict(peak_id='P1',sample='样品',rt_min=5.5)]
        spectrum=np.zeros(1201);spectrum[71]=60;spectrum[89]=100
        self.queries={'P1_BG':spectrum.copy(),'P1_RAW':spectrum.copy()}
        self.hits={'P1_BG':[hit()], 'P1_RAW':[hit()]}
        self.index=[]
        self.evidence=[reference()]

    def run_review(self,**options):
        return summarize(self.peaks,self.queries,link_spectra(self.peaks,self.queries),self.hits,
                         [(14,5),(15,6)],self.index,[],self.evidence,column='HP-5ms',**options)

    def test_missing_adams_retained_with_ions_and_citation(self):
        result=self.run_review();row=result['candidates'][0]
        self.assertEqual(len(result['candidates']),1)
        self.assertIn('精油库里面搜索不到',row['coverage'])
        self.assertEqual(row['name'],'Synthetic candidate (E)-')
        self.assertEqual(row['grade'],'较强暂定候选')
        self.assertEqual(row['sample_ions'][0],dict(mz=89,percent=100.0))
        self.assertIsNone(row['external_evidence'][0]['measured_ions'][2]['percent'])
        self.assertEqual(row['mode'],'BG')

    def test_strict_thresholds_and_low_mf_not_high_confidence(self):
        for mf,rmf,ri,count in [(800,890,'1451',1),(801,800,'1451',0),(801,801,'1470',0)]:
            self.hits={'P1_BG':[hit(mf=mf,rmf=rmf,ri=ri)]}
            rows=self.run_review()['candidates']
            self.assertEqual(len(rows),count)
            if mf==800:self.assertIn('MF不足',rows[0]['grade'])

    def test_raw_fallback_and_same_peak_competitors(self):
        self.hits={'P1_BG':[hit(cas='638-59-5')], 'P1_RAW':[hit(),hit(rank=2)]}
        rows=self.run_review()['candidates']
        self.assertEqual(len(rows),2)
        raw=next(r for r in rows if r['cas']=='6221-98-3')
        self.assertEqual(raw['mode'],'RAW')
        self.assertIn('仅RAW',raw['grade'])
        self.assertEqual(len(raw['raw_bg_hits']['RAW']),2)
        self.assertTrue(any('不能计为' in s for s in raw['caveats']))

    def test_index_presence_outside_window_never_called_absent(self):
        self.index=[dict(name='Present but incompatible',cas='6221-98-3',ai=1500)]
        row=self.run_review()['candidates'][0]
        self.assertEqual(row['coverage'],'精油库检出对应CAS')
        self.assertIn('AI不一致',row['grade'])

    def test_no_extrapolation_and_unverified_reference(self):
        self.peaks[0]['rt_min']=4.9
        row=self.run_review()['candidates'][0]
        self.assertIsNone(row['ri']);self.assertIsNone(row['delta_ri'])
        self.assertIn('无法内插',row['grade'])
        self.peaks[0]['rt_min']=5.5
        self.hits={'P1_BG':[hit(ri='0')]}
        self.assertIn('参考RI缺失',self.run_review()['candidates'][0]['grade'])

    def test_external_conflicts_and_incompatible_modes_not_hidden(self):
        self.evidence=[reference(),reference(ri=1500,id='disagreement')]
        row=self.run_review()['candidates'][0]
        self.assertEqual(len(row['external_evidence']),2)
        self.assertIn('冲突',row['grade'])
        self.evidence=[reference(column='DB-Wax'),reference(mode='isothermal',id='iso')]
        row=self.run_review()['candidates'][0]
        self.assertIn('条件未验证',row['grade'])
        self.assertTrue(all(not e['comparable'] for e in row['external_evidence']))

    def test_blank_exclusion_and_comparable_processing(self):
        self.peaks.append(dict(peak_id='B1',sample='CK',rt_min=5.51))
        self.queries['B1_BG']=self.queries['P1_BG']/2
        self.hits['B1_BG']=[hit()]
        result=self.run_review(blank_samples=['CK'])
        self.assertEqual(result['counts']['candidate_rows'],1)
        best=result['candidates'][0]['blank']['best']
        self.assertEqual(best['blank_peak'],'B1');self.assertAlmostEqual(best['signal_ratio'],2)
        with self.assertRaisesRegex(ValueError,'空白样品名'):
            self.run_review(blank_samples=['typo'])
        del self.queries['B1_BG'];self.queries['B1_RAW']=self.queries['P1_RAW']/2
        self.hits={'P1_BG':[hit()]}
        del self.queries['P1_RAW']
        self.assertIsNone(self.run_review(blank_samples=['CK'])['candidates'][0]['blank']['best'])


class ImportTest(unittest.TestCase):
    def test_decimal_nist_ri_preserved(self):
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'report.txt'
            path.write_text('Unknown: P1\nHit 1 : <<Test>>;<<C18H36O2>>; MF: 870; RMF: 890; Prob: 10; CAS:6221-98-3; Mw: 284; Lib: <<synthetic>>; Id: 1; RI: 1451.75.\n',encoding='utf-8')
            self.assertEqual(parse_report(path)['P1'][0]['library_ri'],'1451.75')

    def test_review_then_summarize_retains_blank_and_column(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            def write(name,text):
                (root/name).write_text(text,encoding='utf-8');return str(root/name)
            peaks=write('peaks.csv','peak_id,sample,rt_min\nP1,样品,5.5\nB1,CK,5.51\n')
            msp=write('spectra.msp','Name: P1_BG\nNum Peaks: 2\n71 60;89 100\nName: B1_BG\nNum Peaks: 2\n71 30;89 50\n')
            anchors=write('alkanes.csv','carbon,rt_min\n14,5\n15,6\n')
            report=write('nist.txt','Unknown: P1_BG\nHit 1 : <<Synthetic (E)-candidate>>;<<C18H36O2>>; MF: 870; RMF: 890; Prob: 10; CAS:6221-98-3; Mw: 284; Lib: <<synthetic>>; Id: 1; RI: 1451.\nUnknown: B1_BG\n')
            evidence=write('evidence.json',json.dumps([reference()]))
            self.assertEqual(main(['review','--peaks',peaks,'--spectra',msp,'--alkanes',anchors,
                                  '--nist-report',report,'--evidence-json',evidence,'--no-bundled-evidence',
                                  '--column','HP-5MS','--blank-sample','CK','--out',str(root/'first')]),0)
            self.assertEqual(main(['summarize','--run-dir',str(root/'first'),
                                  '--out',str(root/'second')]),0)
            first=json.loads((root/'first/candidates.json').read_text(encoding='utf-8'))
            second=json.loads((root/'second/candidates.json').read_text(encoding='utf-8'))
            self.assertEqual(first['candidates'],second['candidates'])
            with ZipFile(root/'first/report.xlsx') as archive:
                pairs=archive.read('xl/worksheets/sheet6.xml').decode()
                self.assertEqual(pairs.count('<row'),1)  # Header only: CK is not a biological replicate.
            for suffix in ('xlsx','csv','html'):
                self.assertTrue((root/'second'/('候选汇总.'+suffix)).is_file())
            self.assertEqual(second['settings']['blank_samples'],['CK'])

    def test_evidence_validation_and_bundle(self):
        self.assertGreaterEqual(len(load_evidence()),13)
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'bad.json'
            record=reference(ri=float('nan'));path.write_text(json.dumps([record]),encoding='utf-8')
            with self.assertRaisesRegex(ValueError,'有限正数'):load_evidence(path,False)

    def test_webbook_conditions_and_lookup_failure(self):
        text='''<p>CAS Registry Number: 6221-98-3</p>
<h3>Van Den Dool and Kratz RI, non-polar column, temperature ramp</h3>
<table><tr><th>Column type</th><th>Active phase</th><th>I</th><th>Reference</th><th>Comment</th></tr>
<tr><td>Capillary</td><td>DB-5</td><td>1984.</td><td>Synthetic author, 2000</td><td>3 K/min</td></tr></table>'''
        records=parse_webbook(text,'6221-98-3','https://example.org/test')
        self.assertEqual(records[0]['ri'],1984)
        self.assertEqual(records[0]['temperature_mode'],'programmed')
        with self.assertRaisesRegex(ValueError,'CAS'):parse_webbook(text,'638-59-5','https://example.org/test')
        with tempfile.TemporaryDirectory() as tmp, patch('gcms_ri.webbook.urlopen',side_effect=HTTPError('https://example.org',403,'blocked',None,None)):
            result=lookup(['6221-98-3'],Path(tmp))
            self.assertIn('失败',result[0]['status'])
            self.assertEqual(json.loads((Path(tmp)/'evidence.json').read_text(encoding='utf-8')),[])


if __name__=='__main__':unittest.main()
