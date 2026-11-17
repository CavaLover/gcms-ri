# GC–MS 保留指数与碎片复核

这个命令行工具把 Agilent 样品读取、NIST 检索、实测保留指数计算、Adams AI 与离子碎片复核汇总到同一份 Excel。完整流程需要 Windows、本机 Agilent 读取组件和 NIST MS Search。已有峰表和质谱的人可以从中间步骤开始。

## 安装与试用

```powershell
python -m pip install "gcms-ri[adams] @ https://github.com/CavaLover/gcms-ri/releases/download/v0.1.0/gcms_ri-0.1.0-py3-none-any.whl"
gcms-ri demo --out demo-output
```

示例使用合成数据，输出 `demo-output/review/report.xlsx`。

## 分析样品

```powershell
gcms-ri doctor --agilent-dir "C:\Lab\Agilent\bin" --nist-dir "C:\NIST17\MSSEARCH"
gcms-ri run --sample "C:\Lab\样品1.D" --sample "C:\Lab\样品2.D" `
  --agilent-dir "C:\Lab\Agilent\bin" --nist-dir "C:\NIST17\MSSEARCH" `
  --alkanes "C:\Lab\alkanes.csv" --adams "C:\Lab\Adams-2017.pdf" `
  --nist-library-version NIST17 --out "C:\Lab\分析001"
```

- `--agilent-dir` 指向包含 `MSDChemDataFile.dll` 及其依赖的安装目录。默认 x64，使用 32 位组件时加 `--agilent-arch x86`。
- `alkanes.csv` 填写已经核实的烷烃标定点，两列为 `carbon,rt_min`。碳数连续，RT 单位为分钟。可先用 `gcms-ri standard --help` 导出标准品峰表和 MSP，再人工核实锚点。
- `--adams` 指向本机 Adams 2017 ed. 4.1 PDF。Adams 候选按 **实测 RI ±20** 查附录 III 的 **AI**。CAS 只在之后与 NIST 候选核对。
- 平行样用 RT 和谱图匹配。碎片按各谱基峰归一。Excel 中的 `RI candidates`、`Fragment ions` 和 `Diagnostic EIC` 可按 `reference_id` 关联。不同采集质量范围可用 `--mass-min`、`--mass-max` 指定。

NIST 检索继承本机设置。如果希望谱库结果与 RI 证据独立，请关闭 NIST 的 **Use RI in Search**，并保留原始命中报告。每次只提交一个任务。

## 用已有结果复核

```powershell
gcms-ri review --peaks 分析001/peaks.csv --spectra 分析001/queries.msp `
  --nist-report 分析001/nist/results.txt --alkanes alkanes.csv `
  --adams Adams-2017.pdf --out 分析002
```

自己导出时，峰表至少包含 `peak_id,sample,rt_min`，MSP 的 `Name` 使用峰号或加 `_RAW`、`_BG` 后缀。已有参考谱可用 `--reference-json` 导入，格式见英文 README。表中名称是暂定候选，完整流程没有做解卷积或峰面积定量。
