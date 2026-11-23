import json
from pathlib import Path
import tempfile
import unittest
from xml.etree import ElementTree
from zipfile import ZipFile

import numpy as np

from gcms_ri.cli import main, parser
from gcms_ri.spectra import MASS_COUNT, write_msp, write_peak_csv


class StudyTest(unittest.TestCase):
    def test_requires_user_anchors(self):
        with self.assertRaises(SystemExit) as exc:
            parser().parse_args(["study","--design","test.json","--out","output"])
        self.assertEqual(exc.exception.code,2)

    def test_study_blank_recheck_homolog_adams_and_workbook(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            run = root / "合成产物"
            (run / "scans").mkdir(parents=True)
            (run / "nist").mkdir()
            peaks = [
                dict(peak_id="A1",sample="幼虫一",rt_min=5.5),
                dict(peak_id="A2",sample="幼虫一",rt_min=8.03),
                dict(peak_id="A3",sample="幼虫一",rt_min=10.03),
                dict(peak_id="N1",sample="幼虫二",rt_min=8.04),
                dict(peak_id="B1",sample="历史CK",rt_min=4.5),
            ]

            def spectrum(mz, intensity):
                values = np.zeros(MASS_COUNT)
                values[mz] = intensity
                return values

            alkane = {57:100,71:75,85:48,99:12,408:10}
            background = {149:100,84:60,131:37}
            spectra = dict(A1_BG=spectrum(list(background),list(background.values())),
                           A2_BG=spectrum(list(alkane),list(alkane.values())),
                           A3_BG=spectrum([57,71,85,99],[100,75,48,12]),
                           N1_BG=spectrum(list(alkane),list(alkane.values())),
                           B1_BG=spectrum([110,155,171],[100,75,50]))
            write_peak_csv(run/"peaks.csv",peaks)
            write_msp(run/"queries.msp",spectra)
            labels = [("a","幼虫一"),("n","幼虫二"),("b","历史CK")]
            (run/"run.json").write_text(json.dumps(
                {"sample_sources":[dict(label=label,name=name) for label,name in labels]}),encoding="utf-8")
            for label,signals in [
                ("a",[(5.50,background,3),(8.03,alkane,3),(10.03,{57:100,71:75,85:48,99:12},3)]),
                ("n",[(8.04,alkane,2)]),
                ("b",[(5.50,background,1),(4.5,{110:100,155:75,171:50},2)]),
            ]:
                rows = []
                for t in np.arange(4,10.5,.004):
                    values = {53:1}
                    for rt,ions,factor in signals:
                        for mz,intensity in ions.items():
                            values[mz] = values.get(mz,0)+intensity*factor*np.exp(-((t-rt)/.014)**2)
                    rows.append(json.dumps(dict(rt=float(t),mz=list(values),intensity=list(values.values()))))
                (run/"scans"/f"{label}.scans.jsonl").write_text("\n".join(rows),encoding="utf-8")
            wrong = ("Hit 1 : <<Tetratriacontane>>;<<C34H70>>; MF: 877; RMF: 890; Prob: 60; "
                     "CAS:14167-59-0; Mw: 478; Lib: <<synthetic>>; Id: 1; RI: 3400.\n")
            (run/"nist"/"results.txt").write_text("".join(
                f"Unknown: {key}\n"+(wrong if key=="A2_BG" else "") for key in spectra),encoding="utf-8")
            anchors = root/"使用者标样.csv"
            anchors.write_text("carbon,rt_min\n28,6\n29,8\n30,10\n",encoding="utf-8")
            design = root/"设计.json"
            design.write_text(json.dumps({"samples":[
                dict(run_dir=run.name,sample=name,group=group,blanks=[dict(sample="历史CK",scope="historical")])
                for name,group in (("幼虫一","第一组"),("幼虫二","第二组"))
            ]},ensure_ascii=False),encoding="utf-8")
            index = root/"adams.json"
            index.write_text(json.dumps([dict(cas="630-03-5",ai=2900,name="Nonacosane",pdf_page=1)]),encoding="utf-8")
            output = root/"复核"
            self.assertEqual(main(["study","--design",str(design),"--alkanes",str(anchors),
                                   "--adams-index",str(index),"--out",str(output)]),0)
            data = json.loads((output/"study.json").read_text(encoding="utf-8"))
            by_peak = {p["peak_id"]:p for p in data["peaks"]}
            self.assertNotIn("B1",by_peak)
            self.assertTrue(by_peak["A1"]["excluded"])
            match = by_peak["A1"]["blank_matches"][0]
            self.assertIsNone(match["peak_id"])  # Found in CK scans, despite absence from its peak table.
            self.assertAlmostEqual(match["sample_blank_ratio"],3,places=1)
            self.assertEqual(match["scope"],"historical")
            self.assertEqual(by_peak["A2"]["candidate"],"正构烷烃 n-C29")
            self.assertEqual(by_peak["A2"]["adams"]["ai"],2900)
            self.assertEqual(by_peak["A2"]["nist"]["BG"][0]["name"],"Tetratriacontane")
            self.assertFalse(by_peak["A3"]["alkane"]["supported"])
            self.assertIsNone(by_peak["A3"]["ri"])
            self.assertEqual([(p["peak_a"],p["peak_b"]) for p in data["pairs"]],[("S1:A2","S2:N1")])
            self.assertEqual([r["carbon"] for r in data["calibration"]["anchors"]],[28,29,30])
            with ZipFile(output/"分组复核.xlsx") as z:
                self.assertIsNone(z.testzip())
                styles = z.read("xl/styles.xml")
                self.assertIn(b"FFBDDDF8",styles)
                namespace = {"s":"http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
                sheet = ElementTree.fromstring(z.read("xl/worksheets/sheet1.xml"))
                rows = sheet.findall(".//s:row",namespace)
                self.assertEqual(len(rows),6)  # Header, three user anchors, two retained experimental peaks.
                self.assertEqual(len(rows[1].findall("s:c[@s]",namespace)),10)  # Full blue anchor row.


if __name__ == "__main__":
    unittest.main()
