# GC–MS 保留指数与候选复核 CLI

把 Agilent 数据读取、NIST 检索、实验 RI、参考谱碎片、CK 和文献证据汇总为可追溯的中文报告。**Adams 精油库未检出对应条目，不再直接排除 NIST 候选。** 程序保留完整检索名、具体离子、参考 RI 条件和未解决的异构体竞争关系。

当前版本 **v0.3.0**。需要 Python 3.11 或更新版本。原生 Agilent/NIST 流程需要 Windows 和用户自己的本地厂商软件。峰表复核、分组复核、候选汇总和合成演示可以跨平台运行。

## 安装与演示

```powershell
python -m pip install --upgrade "gcms-ri[adams] @ https://github.com/CavaLover/gcms-ri/releases/download/v0.3.0/gcms_ri-0.3.0-py3-none-any.whl"
python -m gcms_ri --version
python -m gcms_ri demo --out demo-output
```

也可使用 `gcms-ri` 命令。如果 Windows 提示找不到命令，使用上述 `python -m gcms_ri`，无需修改 PATH。演示全部采用合成输入和合成文献记录，数值不可用于科研引用。

## 分析原始数据与 CK

```powershell
python -m gcms_ri doctor --agilent-dir "C:\Lab\Agilent\bin" --nist-dir "C:\NIST17\MSSEARCH"
python -m gcms_ri run --sample "C:\Lab\样品1.D" --sample "C:\Lab\样品2.D" `
  --blank "C:\Lab\CK.D" --column HP-5ms `
  --agilent-dir "C:\Lab\Agilent\bin" --nist-dir "C:\NIST17\MSSEARCH" `
  --alkanes "C:\Lab\alkanes.csv" --adams "C:\Lab\Adams-2017.pdf" `
  --nist-library-version NIST17 --out "C:\Lab\分析001"
```

- `--blank` 指定 CK 的 `.D`，可以重复。CK 被读取并单独核查，不参与生物学重复配对。
- `--sample-dir` 只适用于直接子目录均为待分析样品的目录。不要把含烷烃标准、其他实验数据的整个目录直接传入；逐项使用 `--sample` 更明确。
- `--agilent-dir` 包含 `MSDChemDataFile.dll` 及依赖；默认 x64，32 位组件使用 `--agilent-arch x86`。读取器使用 Windows .NET Framework 编译器，只读取原始 `data.ms` 的副本。
- `--column` 填实际色谱柱。未填写时，外部 RI 不会被标成“条件已验证”。固定相相近也不代表升温程序完全相同。
- `--adams` 支持用户本地的 Adams 2017 第 4.1 版 PDF。也可用 `--reference-json` 导入本地参考谱；两者均可省略，此时不会声称已经核查 Adams 覆盖情况。

烷烃表须来自相同 GC 条件下已经人工核实的标准：碳数连续、RT 递增、单位为分钟。示例格式如下，数值仅演示：

```csv
carbon,rt_min
14,5.000
15,6.000
```

可先执行 `python -m gcms_ri standard --standard 标准.D --agilent-dir "C:\Lab\Agilent\bin" --out 标准复核`，导出候选时间与 MSP，再核实碳数。程序不自动把某个峰指定为某个正构烷烃。

NIST 检索沿用本机配置。为使谱库检索与 RI 筛选相互独立，请关闭 **Use RI in Search**，并导出含 MF、RMF、CAS、库编号和 RI 的命中报告。程序保留设置快照。若只导出每谱前两项，只能复核这两项；建议导出更多命中供竞争候选审查。NIST 共用后台文件，每次只运行一个检索任务；超时后须先确认 NIST 空闲，再处理保留的锁文件。

## 分组复核

`study` 读取已有的 `run` 或 `review` 目录。每份实验样品在设计文件中填写分组和对应 CK。保留 `run` 目录中的 `run.json` 和 `scans` 可让程序复核提取离子色谱（EIC）。相对 `run_dir` 以设计文件所在目录为基准。

```json
{
  "samples": [
    {
      "run_dir": "本批分析",
      "sample": "样品A",
      "group": "处理组",
      "blanks": [{"sample": "CK", "scope": "same_batch"}]
    },
    {
      "run_dir": "另一批分析",
      "sample": "样品B",
      "group": "比较组",
      "blanks": [{"run_dir": "本批分析", "sample": "CK", "scope": "historical"}]
    }
  ]
}
```

```powershell
python -m gcms_ri study --design "C:\Lab\分组.json" --alkanes "C:\Lab\本机实测标样.csv" --out "C:\Lab\分组复核"
```

- **`--alkanes` 是必填项。** 提供自己核实过的 `carbon,rt_min` 校准 CSV，或先用 `standard --standard` 读取自己的原始 Agilent `.D` 再整理校准文件。蓝色行使用该 CSV 的保留时间和碳数。多个批次共用校准前，需要确认采集方法和标准峰可比。跨条件的数据应各自使用对应校准分开复核。
- 每份样品可指定多个 CK，`scope` 可填 `same_batch`、`historical` 或 `unspecified`。跨批次 CK 用于排查代表性背景。默认比较保留时间差 ≤0.10 min 和同处理模式 EI 余弦 ≥0.80。提供原始扫描时，再通过特征离子提取补查 CK 中未自动检出的峰。
- 样品和 CK 的峰形相符时输出局部基线校正后的特征离子面积比。面积用于衡量空白贡献，结果不做 CK 直接扣除或浓度换算。同峰及典型硅氧烷、邻苯二甲酸酯类背景候选移出分组主表，全部证据留在 `CK审计`。
- 正构烷烃结合实测标样位置、烷烃碎片及分子离子同峰暂定归属。可在确认漂移后显式填写每个样品的 `rt_offset_min`，供匹配窗口使用。样品实验 RI 始终用实测 TR 内插。标定范围外不外推。
- Adams 对应条目按候选 CAS 单列，来源为运行目录缓存的 Adams 索引或显式传入的 `--adams-index`。输出保留 NIST 原始命中，方便复核同系物排序冲突。通过标样和分子离子得到的归属仍是暂定。

输出为 `分组复核.xlsx` 和 `study.json`。分组主表使用相同表头，蓝色标样按 TR 穿插。另有样品与 CK、分组烷烃对照、样品间峰配对和 CK 审计页。原有 `run`、`review` 与 `summarize` 的输出仍适用于原工作流。

## 复用旧结果

已经有 `peaks.csv`、`queries.msp`、`alkanes.csv` 和 `nist/results.txt` 时，无需重新进行 NIST 检索：

```powershell
python -m gcms_ri summarize --run-dir "C:\Lab\分析001" `
  --blank-sample CK --column HP-5ms --out "C:\Lab\候选复核002"
```

`--blank-sample` 必须与峰表 `sample` 列的完整名称一致，可重复。v0.2.0 输出会缓存完整参考索引和窗口内参考谱，汇总时自动复用。没有另给文献参数时也复用 `evidence-used.json`，保留人工核查记录；重新指定 `--evidence-json` 可更新证据。旧版没有缓存时，程序尝试读取 `run.json` 中仍存在的本地参考 PDF。

也可显式提供完整索引，避免再次解析 PDF：

```powershell
python -m gcms_ri summarize --run-dir "C:\Lab\旧分析" `
  --reference-index "C:\Lab\完整Adams索引.json" --index-kind Adams `
  --blank-sample CK --column HP-5ms --out "C:\Lab\新汇总"
```

索引格式为 `[ {"name":"名称","cas":"CAS号","ai":1500,"pdf_page":100} ]`。必须提供完整索引；只有 RI 窗口子集不能用于判断“精油库未检出”。仅有索引时，报告不会虚构对应参考谱的离子强度。自定义参考表默认 `--index-kind local`，不冒充 Adams。

已有导出文件，也可重做完整的 RI 窗口和碎片比较：

```powershell
python -m gcms_ri review --peaks 分析001/peaks.csv --spectra 分析001/queries.msp `
  --nist-report 分析001/nist/results.txt --alkanes 分析001/alkanes.csv `
  --adams Adams-2017.pdf --column HP-5ms --blank-sample CK --out 分析003
```

峰表至少含 `peak_id,sample,rt_min`；MSP 的 `Name` 使用唯一 ASCII 峰号，或加 `_RAW`、`_BG` 后缀。参考谱 JSON 示例：

```json
[{"name":"Synthetic reference","cas":"121-44-8","ai":1508,"mass_min":38,"mass_max":180,"ions":{"43":66,"55":75,"69":100,"83":73}}]
```

## 筛选规则

1. 同峰同 CAS 合并，优先采用满足纳入条件的 BG 命中，其次采用 RAW，保留两种谱的分数与其他命中。
2. 默认 `RMF > 800`，且 `|实验 RI − 初筛参考 RI| < 20`。先核对相同 CAS 的 Adams/本地参考 AI；未匹配时保留 NIST 库 RI 作为初筛线索，不因 Adams 缺条目而删除。
3. “较强暂定候选”还要求 `MF > 800`、BG 支持，并有本地参考 AI 或相近固定相的程序升温 RI 支持。外部可比 RI 有冲突时降级；全部来源都保留，不能挑一个最接近的值掩盖冲突。
4. 高 MF/RMF 但实验 RI 无法内插或参考 RI 缺失的记录，保留为“待复核”。`--keep-peak 峰号` 可额外保留该峰 RMF 过阈值的命中，RI 不一致会明确标记。
5. CK 以 RT 差 ≤0.10 min、同处理模式的谱余弦 ≥0.80 匹配。未匹配不等于空白不存在；峰顶信号比不是浓度比，也没有直接进行 CK 扣除。

阈值可用 `--mf-min`、`--rmf-min`、`--ri-window` 调整。原详细报告的独立 RI 窗口表采用包含边界的 `≤`；候选汇总采用上述严格 `<`，两表用途不同。

**RMF 不是鉴定准确率。** 程序不做共洗脱解卷积、不报告峰面积或浓度。主要离子是指定质量范围内最强 8 个信号，不自动等于诊断离子。实验 RI 只在烷烃锚点之间线性插值，范围外不外推；非手性 RI 与 EI 不能单独确定对映体。所有结果仍为暂定鉴定。

## 文献与外部 RI

附带 13 条已核查的公开 RI/文献记录，包括丁酸十四烷基酯、乙酸十四烷基酯、长链醛及相关竞争候选。只包含少量书目信息、RI 和文献列出的离子编号，**不包含 NIST/Adams 商业谱库或原始实验数据**。

附带记录区分已读取原文/补充材料、经 NIST 数据库追溯和数据库保留数据。这不是穷尽性文献搜索；未命中记录会明确显示“本轮未检出已核实的对应文献记录”。

按 CAS 联网查询 NIST Chemistry WebBook 的公开 RI 表与文献线索：

```powershell
python -m gcms_ri literature --cas 6221-98-3 --cas 638-59-5 --out 文献检索001
python -m gcms_ri summarize --run-dir 分析001 --column HP-5ms `
  --evidence-json 文献检索001/evidence.json --out 文献复核004
```

联网只发送 CAS，不上传样品或谱图。网页快照、校验值及检索日志写入输出目录。网络拒绝或页面格式变化会记录“检索失败”，不能解释为没有文献；部分失败返回非零退出码，并保留成功结果。自动提取的记录不代表读过原论文全文。

可手工准备 `--evidence-json`，见 [字段与证据说明](docs/evidence-format.md)。`--no-bundled-evidence` 可关闭附带记录。人工新增文献须自行核读，并如实填写核查范围。

## 输出

- `候选汇总.xlsx`：中文候选总表、碎片与命中、文献索引、筛选说明，可筛选。
- `候选汇总.html`：可按样品、名称、CAS、状态搜索。
- `候选汇总.csv`：中文字段，UTF-8 BOM。
- `candidates.json`：完整数值、逐候选判定、CK、竞争命中和外部证据，供继续分析。
- `report.xlsx`：`run/review` 的原有详细技术报告，含所有峰、NIST、独立 RI 窗口、碎片、重复匹配及 EIC；保留旧表名兼容既有脚本。
- `run.json`、`reference-index.json`、`reference-spectra.json`、`evidence-used.json`：参数、来源与本地缓存。缓存谱图可能来自用户自备的商业参考资料，应按其授权使用。

输出目录必须为新目录，原始数据不覆盖。

## 开发与测试

```powershell
python -m pip install -e ".[adams,test]"
python -m unittest discover -s tests -v
python -m pip install build
python -m build
```

GitHub Actions 在 Windows/Python 3.13 与 Linux/Python 3.11 上测试与打包。仓库采用 MIT 许可，仅分发工具代码、公开书目信息和合成测试。Agilent 组件、NIST 软件/谱库和 Adams PDF 须由用户自行合法获取。
