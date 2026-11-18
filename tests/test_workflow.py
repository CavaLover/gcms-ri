import json
from pathlib import Path
import tempfile
import unittest
from xml.etree import ElementTree
from zipfile import ZipFile

from gcms_ri.cli import main
from gcms_ri.spectra import interpolate, load_anchors, parse_msp


class WorkflowTest(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.root = Path(tmp.name)

    def write(self, name, contents):
        path = self.root / name
        path.write_text(contents, encoding="utf-8")
        return str(path)

    def test_demo_from_command_to_xlsx(self):
        output = self.root / "合成 示例"
        self.assertEqual(main(["demo", "--out", str(output)]), 0)
        with ZipFile(output / "review" / "report.xlsx") as archive:
            self.assertIsNone(archive.testzip())
            main_xml = ElementTree.fromstring(archive.read("xl/workbook.xml"))
            sheet_ns = {"s": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
            names = [sheet.attrib["name"] for sheet in main_xml.findall(".//s:sheet", sheet_ns)]
            self.assertEqual(names[:5], ["Summary", "Peaks", "NIST", "RI candidates", "Fragment ions"])
            candidates = archive.read("xl/worksheets/sheet4.xml").decode("utf-8")
            self.assertIn("Synthetic reference A", candidates)
            self.assertNotIn("Synthetic reference outside the RI window", candidates)
            self.assertEqual(candidates.count("<row"), 5)

    def test_nist_cas_cannot_add_or_remove_ri_candidates(self):
        output = self.root / "review"
        peaks = self.write("peaks.csv", "peak_id,sample,rt_min\nP1,S,5.5\nP2,S,7\n")
        msp = self.write("spectra.msp", """Name: P1_BG
Num Peaks: 2
43 100; 57 90
Name: P2
Num Peaks: 2
43 100
57 90
""")
        alkanes = self.write("alkanes.csv", "carbon,rt_min\n14,5\n15,6\n")
        refs = self.write("references.json", json.dumps([
            {"name": "Exactly on window boundary", "ai": 1470, "ions": {43: 100, 57: 90}},
            {"name": "Independent NIST match", "ai": 1454, "cas": "121-44-8", "ions": {43: 100, 57: 90}},
            {"name": "Outside window despite NIST", "ai": 1471, "cas": "112-92-5", "ions": {43: 100, 57: 90}},
        ]))
        hits = self.write("nist.txt", """Unknown: P1_BG
Hit 1  : <<Mock hit A>>;<<C18H38O>>; MF: 911; RMF: 933; Prob: 10.25; CAS:112-92-5; Mw: 270; Lib: <<demo>>; Id: 1; RI: 1471.
Hit 2  : <<Mock hit B>>;<<C6H15N>>; MF: 800; RMF: 900; Prob: 5.00; CAS:121-44-8; Mw: 101; Lib: <<demo>>; Id: 2; RI: 1454.
Unknown: P2
""")
        self.assertEqual(main(["review", "--peaks", peaks, "--spectra", msp, "--alkanes", alkanes,
                               "--reference-json", refs, "--nist-report", hits, "--out", str(output)]), 0)
        with ZipFile(output / "report.xlsx") as package:
            data = package.read("xl/worksheets/sheet4.xml").decode()
        self.assertIn("Exactly on window boundary", data)
        self.assertIn("BG #2", data)
        self.assertNotIn("Outside window despite NIST", data)
        self.assertEqual(data.count("<row"), 3)

    def test_anchor_boundary_and_bad_msp_count(self):
        path = Path(self.write("alkanes.csv", "carbon,rt_min\n14,5\n15,6\n"))
        anchors = load_anchors(path)
        self.assertEqual(interpolate(5.5, anchors), 1450)
        self.assertEqual(interpolate(6, anchors), 1500)
        self.assertIsNone(interpolate(6.01, anchors))
        with self.assertRaisesRegex(ValueError, "Num Peaks"):
            parse_msp(Path(self.write("bad.msp", "Name: broken\nNum Peaks: 2\n43 -50\n")))


class AdamsLayoutTest(unittest.TestCase):
    def test_vector_fragment_extraction_from_synthetic_pdf(self):
        from reportlab.pdfgen import canvas

        from gcms_ri.adams import pdf_entries

        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "synthetic.pdf"
            page = canvas.Canvas(str(path), pagesize=(612, 792))
            page.setFont("Helvetica", 10)
            page.drawString(36, 700, "RT: 35.64 AI: 1671 KI: 1672 Synthetic alcohol")
            page.drawString(36, 686, "CAS#: 112-72-1")
            origin_x, baseline, y_height = 54, 540, 90
            scale = 1.80
            path_shape = page.beginPath()
            path_shape.moveTo(origin_x, baseline + y_height)
            path_shape.lineTo(origin_x, baseline)
            path_shape.lineTo(origin_x + scale * (220 - 38), baseline)
            page.drawPath(path_shape)
            for mass in range(40, 221, 20):
                x = origin_x + (mass - 38) * scale
                page.line(x, baseline, x, baseline - 3)
                page.drawCentredString(x, baseline - 12, str(mass))
            for mass, abundance in ((41, 100), (55, 85), (83, 72), (196, 2)):
                x = origin_x + (mass - 38) * scale
                high = baseline + y_height * abundance / 100
                page.line(x, baseline, x, high)
                page.drawCentredString(x, high + 3, str(mass))
            page.showPage()
            page.setFont('Helvetica',10)
            page.drawString(36,700,'RT: 45.64 AI: 2000 KI: 2000 Synthetic outside window')
            page.drawString(36,686,'CAS#: 638-66-4')
            page.save()
            index=[]
            references, count = pdf_entries(path, [1677.3], 20, index)
            self.assertEqual(count, 2)
            self.assertEqual(len(index),2)
            self.assertEqual(len(references),1)
            self.assertEqual(index[1]['cas'],'638-66-4')
            self.assertEqual(references[0]["ai"], 1671)
            self.assertEqual(references[0]["spectrum_status"], "Compared")
            self.assertAlmostEqual(references[0]["ions"][196], 2, places=2)


if __name__ == "__main__":
    unittest.main()
