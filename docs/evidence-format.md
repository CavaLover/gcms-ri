# 文献及外部 RI 格式

`--evidence-json` 接收 JSON 数组，使用 UTF-8 编码。示例为真实公开文献的最小事实记录：

```json
[
  {
    "id": "my-review-oliveira-2015-s2",
    "cas": "6221-98-3",
    "name": "Tetradecyl butyrate",
    "ri": 1984,
    "column": "DB-5ms",
    "temperature_mode": "programmed",
    "citation": "Oliveira et al. (2015), The origin and evolution of queen and fertility signals in Corbiculate bees, Table S2",
    "url": "https://doi.org/10.1186/s12862-015-0509-8",
    "source_url": "https://doi.org/10.1186/s12862-015-0509-8",
    "verification": "已读取原始补充材料",
    "note": "补充表实验RI1984、库RI1985；仍需排除异丁酸酯竞争候选。",
    "ions": [71, 89, 168, 196, 241]
  }
]
```

`id`、`cas`、`citation`、`url`、`column`、`temperature_mode`、`verification` 必填且不为空。编号不能与附带记录重复；附带编号以 `curated-` 开头。`ri` 是正数或 `null`；只有名称报道而没有可比 RI 时，使用 `null`。`ions` 是文献列出的整数 m/z，不是虚构的相对强度。工具另列这些 m/z 在样品中的实际强度；采集范围外显示未获得，不填零。

| temperature_mode | 含义 |
|---|---|
| `programmed` | 程序升温 |
| `isothermal` | 恒温 |
| `unknown` | 未核实 |

工具只把已识别的相近固定相族、且均为程序升温的数据用于自动外部 RI 支持判断。当前识别常见的 5% 苯基柱、甲基硅氧烷柱和 Wax 柱；未知柱名不自动归类。不能把预测 RI 写成实测程序升温 RI；预测或条件未核实的数据应填写 `unknown` 并在 `note` 说明。

同一 CAS 可以有多条不同编号的文献。程序保留所有条目、计算每条差值，并显示可比记录中的冲突。数据库 RI 与本地 NIST 库 RI 可能来自同一文献，不应称为两个独立实验验证。

`literature` 子命令从 NIST 公开网页提取条件和引用线索，保存 HTML 及 SHA-256。它不会下载商业参考谱，也不会自动断言论文讨论了蚂蚁或信息素功能。请人工核查 DOI、原始表格、固定相及目标异构体。
