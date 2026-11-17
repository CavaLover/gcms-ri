"""Generate a redistributable, explicitly synthetic GC–MS example."""

import json


def demo(out, execute):
    from .cli import create_output, parser

    create_output(out)
    fixture = out / "synthetic_input"
    fixture.mkdir()
    peaks = fixture / "peaks.csv"
    peaks.write_text("peak_id,sample,rt_min\ns1_P001,Demo A,6.100\ns2_P001,Demo B,6.105\n", encoding="utf-8")
    alkane = fixture / "alkanes.csv"
    alkane.write_text("carbon,rt_min\n14,5.0\n15,6.0\n16,7.0\n", encoding="utf-8")
    spectra = fixture / "spectra.msp"
    spectra.write_text(
        "Name: s1_P001_BG\nNum Peaks: 5\n43 60\n55 80\n69 100\n83 75\n168 3\n\n"
        "Name: s2_P001_BG\nNum Peaks: 5\n43 57\n55 82\n69 100\n83 77\n168 3.5\n", encoding="ascii"
    )
    references = fixture / "references.json"
    references.write_text(json.dumps([
        {"name": "Synthetic reference A", "ai": 1508, "mass_min": 38, "mass_max": 180,
         "ions": {"43": 66, "55": 75, "69": 100, "83": 73, "168": 2.8}},
        {"name": "Synthetic reference B", "ai": 1522, "mass_min": 38, "mass_max": 180,
         "ions": {"43": 100, "59": 80, "71": 35, "106": 10}},
        {"name": "Synthetic reference outside the RI window", "ai": 1550, "mass_min": 38, "mass_max": 180,
         "ions": {"43": 60, "55": 80, "69": 100}},
    ], indent=2), encoding="utf-8")
    args = parser().parse_args([
        "review", "--peaks", str(peaks), "--spectra", str(spectra), "--alkanes", str(alkane),
        "--reference-json", str(references), "--out", str(out / "review"),
    ])
    print("Synthetic inputs; no vendor software or real NIST/Adams data is used.", flush=True)
    execute(args)
