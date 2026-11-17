"""XLSX writer using public, installable dependencies."""

import xlsxwriter


PEAKS = ("peak_id", "sample", "rt_min", "ri", "tic", "primary_query")
HITS = ("query", "rank", "name", "formula", "mf", "rmf", "prob", "cas", "mw", "library", "library_id", "library_ri")
REFERENCES = ("reference_id", "peak_id", "sample", "rt_min", "ri", "name", "cas", "ai", "delta_ri", "pdf_page",
              "plot_type", "spectrum_status", "mass_min", "mass_max", "spectrum", "cosine", "nist_same_cas")
IONS = ("reference_id", "peak_id", "mz", "reference_percent", "measured_bg_percent", "measured_raw_percent")


def book(path, peaks, hits, matches, fragments, pairs, diagnostic, anchors, metadata):
    with xlsxwriter.Workbook(path, {"strings_to_formulas": False, "strings_to_urls": False,
                                   "constant_memory": True}) as workbook:
        header = workbook.add_format({"bold": True, "font_color": "#FFFFFF", "bg_color": "#254D6B",
                                      "valign": "vcenter"})
        numerics = {
            "rt_min": workbook.add_format({"num_format": "0.0000"}),
            "delta_rt_min": workbook.add_format({"num_format": "0.0000"}),
            "ri": workbook.add_format({"num_format": "0.0"}),
            "delta_ri": workbook.add_format({"num_format": "0.0"}),
            "cosine": workbook.add_format({"num_format": "0.000"}),
            "pearson_r": workbook.add_format({"num_format": "0.000"}),
        }
        summary = [
            ("Report", "GC–MS tentative candidate review"),
            ("Peaks", len(peaks)), ("Peaks with interpolated RI", sum(p["ri"] is not None for p in peaks)),
            ("NIST queries", len(hits)), ("Adams/reference matches within RI window", len(matches)),
            ("Compared reference spectra", sum(r["cosine"] is not None for r in matches)),
            ("Replicate pairs", len(pairs)),
            ("Reference selection", "Reference AI is selected using measured RI and the configured window."),
            ("NIST CAS", "Shown after independent RI selection. MF/RMF and Prob are library report values."),
            ("Status", "All assignments are tentative. No deconvolution, peak area or concentration is reported."),
        ]
        for key, value in metadata.items():
            if isinstance(value, (str, int, float)):
                summary.append((key, value))
        _table(workbook, "Summary", ("Field", "Value"), summary, (46, 118), header, numerics)
        _table(workbook, "Peaks", PEAKS, _values(peaks, PEAKS), (19, 28, 15, 12, 15, 22), header, numerics)
        hit_records = [dict(query=query, **hit) for query, records in hits.items() for hit in records]
        _table(workbook, "NIST", HITS, _values(hit_records, HITS),
               (21, 9, 51, 17, 9, 9, 12, 17, 11, 15, 13, 13), header, numerics)
        _table(workbook, "RI candidates", REFERENCES, _values(matches, REFERENCES),
               (16, 19, 27, 15, 12, 54, 17, 11, 13, 12, 14, 44, 13, 13, 13, 13, 29), header, numerics)
        _table(workbook, "Fragment ions", IONS, fragments, (16, 19, 11, 24, 24, 24), header, numerics)
        _table(workbook, "Replicate pairs", ("sample_a", "sample_b", "peak_a", "peak_b", "delta_rt_min", "cosine"),
               pairs, (26, 26, 19, 19, 16, 14), header, numerics)
        _table(workbook, "Diagnostic EIC", ("reference_id", "peak_id", "reference_mz", "target_mz", "pearson_r", "target_apex_shift_sec"),
               diagnostic, (16, 19, 18, 17, 15, 24), header, numerics)
        _table(workbook, "Calibration", ("carbon", "rt_min", "ri"),
               [(n, rt, n * 100) for n, rt in anchors], (14, 16, 14), header, numerics)


def _values(rows, keys):
    return (tuple(record[k] for k in keys) for record in rows)


def _table(workbook, name, columns, rows, widths, header, num_formats):
    sheet = workbook.add_worksheet(name)
    sheet.freeze_panes(1, 1)
    sheet.hide_gridlines(2)
    sheet.set_row(0, 22)
    for i, (column, size) in enumerate(zip(columns, widths, strict=True)):
        sheet.set_column(i, i, size, num_formats.get(column))
    sheet.write_row(0, 0, columns, header)
    row_number = 0
    for row_number, values in enumerate(rows, 1):
        if row_number >= 1_048_576:
            raise ValueError(f"Excel row limit exceeded on {name}. Use a narrower RI window or smaller sample batch.")
        sheet.write_row(row_number, 0, values)
    sheet.autofilter(0, 0, row_number, len(columns) - 1)
