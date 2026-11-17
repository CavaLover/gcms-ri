# GC–MS RI Review

GC–MS RI Review is a Python command-line tool for documenting tentative GC–MS identifications. It can extract apex spectra from Agilent data on Windows, submit MSP to an installed NIST MS Search, interpolate measured retention indices, and compare candidates against local Adams AI and reference EI spectra.

[中文使用说明](README.zh-CN.md)

## Install and try

Requires Python 3.11 or newer. Install the wheel in your own Python environment.

```powershell
python -m pip install "https://github.com/CavaLover/gcms-ri/releases/download/v0.1.0/gcms_ri-0.1.0-py3-none-any.whl"
gcms-ri --version
gcms-ri demo --out demo-output
```

The demo creates synthetic spectra and an XLSX at `demo-output/review/report.xlsx`. It works without NIST, Agilent, or the Adams PDF.

Enable Adams PDF parsing with:

```powershell
python -m pip install "gcms-ri[adams] @ https://github.com/CavaLover/gcms-ri/releases/download/v0.1.0/gcms_ri-0.1.0-py3-none-any.whl"
```

## Analyze Agilent .D on Windows

The native reader needs a local Agilent installation containing `MSDChemDataFile.dll` and its dependencies, plus the Windows .NET Framework C# compiler. Supply the vendor DLL directory to `--agilent-dir`. The `--agilent-arch` switch selects `x64` (default) or `x86`. The NIST workflow was validated with NIST 17 MS Search. Set NIST Search to export hit lines that include MF, RMF, CAS, library ID and RI. For independent evidence, disable **Use RI in Search** in NIST. Each run snapshots the current NIST settings.

```powershell
gcms-ri doctor --agilent-dir "C:\Lab\Agilent\bin" --nist-dir "C:\NIST17\MSSEARCH"
gcms-ri run --sample "C:\Lab\sample-A.D" --sample "C:\Lab\sample-B.D" `
  --alkanes "C:\Lab\alkanes.csv" --agilent-dir "C:\Lab\Agilent\bin" `
  --nist-dir "C:\NIST17\MSSEARCH" --nist-library-version "NIST17" `
  --adams "C:\Lab\Adams-2017.pdf" --out "C:\Lab\review-001"
```

Use `--sample-dir` for a directory containing only sample `.D` folders. The native mode reports peak spectra, not peak areas, and reads copies of `data.ms`. The NIST bridge uses the installed desktop program, its shared background-search files and the local configuration. Run one search at a time and keep other NIST/ChemStation submissions idle. Check NIST is idle before clearing a lock left after a timeout. If the staging path contains non-ASCII characters, give `--nist-temp` an existing writable ASCII-only directory.

Provide a confirmed alkane calibration table with consecutive carbon numbers and retention times in minutes. The two rows below illustrate the format.

```csv
carbon,rt_min
14,5.000
15,6.000
```

To export candidate times and MSP spectra from a standard injection before assigning carbon numbers, use `gcms-ri standard --standard standard.D --agilent-dir "C:\Lab\Agilent\bin" --out standard-check`. Use standards acquired with the same GC method as the samples.

## Import results and rerun review

```powershell
gcms-ri review --peaks review-001/peaks.csv --spectra review-001/queries.msp `
  --nist-report review-001/nist/results.txt --alkanes alkanes.csv `
  --adams Adams-2017.pdf --ri-window 20 --out review-002
```

The peak CSV must have `peak_id,sample,rt_min`. Use unique ASCII peak IDs with no spaces. MSP `Name` must be the peak ID, or the ID suffixed by `_RAW` and/or `_BG`. If you have no NIST results, omit `--nist-report`. You can also generate a report with `gcms-ri nist queries.msp --nist-dir "C:\NIST17\MSSEARCH" --out search-001`.

`--reference-json` can be used in place of `--adams` with a locally prepared EI reference list. Example:

```json
[{"name":"Synthetic reference","ai":1508,"mass_min":38,"mass_max":180,"ions":{"43":66,"55":75,"69":100,"83":73,"168":2.8}}]
```

The optional fields `cas` and `pdf_page` are retained in the output. The PDF parser supports the text and plot layout of Adams 2017, ed. 4.1, Appendix III. When a raster axis cannot be verified, the RI candidate stays in the sheet with a review status. Use the reference document and relevant standards to adjudicate weak or ambiguous spectra.

## Read the result

- **Peaks / Calibration** show peak RT, measured RI, and alkane anchors. RI is interpolated only between available anchors.
- **NIST / RI candidates** show independent library hits and reference candidates. Adams candidates are selected solely by `|measured RI − AI| ≤ 20` by default. The report then annotates any matching NIST CAS numbers.
- **Fragment ions** shows m/z and abundance normalized within each spectrum. Cosine is computed over the common m/z window. Set `--mass-min` and `--mass-max` to your acquisition range (default 38–600).
- **Replicate pairs / Diagnostic EIC** show RT and spectrum matches across samples and targeted EIC checks for the C14/C16/C18 1-alkanol reference entries when native scans are available.
- **run.json** records input hashes, settings and stage parameters. `peaks.csv`, `queries.msp` and the NIST report can be reused for later review.

This workflow supports tentative identification. It performs local background subtraction and does not deconvolve coeluting compounds. Standard/filter choices and detector thresholds can change which peaks are detected.

The repository distributes the MIT-licensed tool and synthetic example only. Obtain your own copy of the third-party components, reference PDF, libraries, and spectrum files.

## Develop

```powershell
python -m pip install -e ".[adams]"
python -m unittest discover -s tests
python -m build
```
