# 附手册天然铁靶样本 Fe-M01 原生输入与读回摘录

本目录的五个原生输入与宿主选项文件取自 Fe-M01 第二轮回答：`geometry.gdml`、`beam.mac`、`physics.txt` 与 `host-options.json` 为 Geant4 侧的几何、源、物理表及宿主选项，`fluka.inp` 为完整 FLUKA 输入。每个文件的 SHA256 与[首轮记录](../records/Fe-M01-r1.json)及[第二轮记录](../records/Fe-M01-r2.json)中同名字段的摘要完全一致，证实首、次轮之间五个原生输入与宿主字段逐字未变。

[读回摘录](readback.json)保留 **`START 1.0` 初始化运行输出**（`Fe-M01-r1` 零事件 Geant4 初始化与 1 个初级 FLUKA 运行）的原文与行号，对应正文 §2 表 2 与 §3.1 误报分析的检查核判定字段：

1. **`incl-physics.tsv`（Geant4 级联团簇开关 P03）**：第 2、3、6、17、18 行记录 `mode=default`、`phase=idle_before_first_BeamOn`、`configuration_verified=1`、`before_cluster=0`、`after_cluster=0`，确认首次束流前 INCL++ 默认团簇算法处于开启状态。
2. **`m1-readback.tsv`（Geant4 产生阈 P24、反冲产生射程 P26、名义温度 P29、中子过程名 P23）**：第 11–13 行 γ/e± 产生射程 `target_{gamma,electron,positron}_cut_mm_bits` 均为 `3fb999999999999a`（IEEE-754 双精度 0.1 mm，对应表 2 P24）；第 14 行 `target_proton_cut_mm_bits` 为 `0000000000000000`（0.0 mm，与 γ/e± 默认 0.1 mm 分开，对应 §3.1 P26）；第 19 行 `material_temperature_K_bits` 为 `4072800000000000`（296.0 K）；第 302–307 行列出全部 6 个中子过程名，均不含 `killer` 或 `cut`。
3. **`input001.out` 与 `fluka.inp`（FLUKA 级联聚结 P03、电磁产生阈 P24、中子输运截止 P23、计分汇总 S05–S07）**：`fluka.inp` 第 35–48 行定义 7 个编码 `124.0` 的 `USRYIELD` 计分器（`FE01`–`FE07`，延续卡各设 4 个极角区、`ixa=3`、格式化单元 `+21`，合计 28 格）；`input001.out` 第 142 行回显 `PHYSICS 1.0 ... COALESCE`，第 827–830 行读回靶材 `NATTARG` 电磁产生阈 `Ae = 0.611000 MeV` 与 `Ap = 0.100000 MeV`（对应表 2 P24），第 1055 行读回中子输运动能下限 `1.000E-14 GeV`（即 $10^{-11}\ \mathrm{MeV}$），第 1245–1248 行确认首个 `USRYIELD` 探测器 `FE01` 的粒子、穿越面、4 个极角区、用户归一 `1.0000E+00` 与采用截面因子 `1.0000E+00 mb`。（28 个编码 `2401.0` 的单格计分器属于开发样本 Fe-M02，见[早期计分案例](../records/early-Fe-scoring.json)。）

## 正式运行原生示例：`formal-g000-f00/`

与上述 `START 1.0` 初始化摘录区分，[formal-g000-f00/](formal-g000-f00/pair-p000-provenance.json) 绑定正式铁靶输运的首对记录 **`Fe-M01/g000-f00`**（对应 `units/Fe-M01.json` 两侧数组下标 `0`）：**Geant4 第 0 批**（`u000`，种子 `70690000`，`10,000` 初级质子）与 **FLUKA 第 0 流**（`u000`，种子 `72690000`，`200,000` 初级质子／总权重 `200,000`）各自统计独立，并非同一批粒子历史，仅在成对记录 `p000` 中成对核查。

- **索引与核查依据绑定**：[重查索引](formal-g000-f00/kernel-index.json)（SHA256 与 `units/Fe-M01.json` 中 `kernel_index_sha256` 完全一致，区分于早期原始目录旧索引）首条指向 `p000/kernel.json`，[p000 来源与哈希摘录](formal-g000-f00/pair-p000-provenance.json)绑定该文件 SHA256、全部 24 个 `runtime.files` 摘要（省略 `kernel.json`/`plan.json` 中的集群本机路径行并标明行号）。
- **FLUKA 第 0 流原件**：[正式输入 input.inp](formal-g000-f00/fluka/input.inp)（在 `fluka.inp` 基础上将 `RANDOMIZ` 与 `START` 分别更新为 `72690000.` 与 `200000.`）、[完整运行输出 input001.out](formal-g000-f00/fluka/input001.out)（第 172、1428 行见种子 `72690000`，第 241、1678–1679 行见 `200000` 初级与总权重 `2.000000E+05`，第 199–238、1243–1299 行见 `FE01`–`FE07` 计分器、角度/动能变量、单位 `21`、用户归一 `1.0000E+00` 与采用截面 `1.0000E+00 mb`）及[完整原生谱 input001_fort.21](formal-g000-f00/fluka/input001_fort.21)（原生格值本身已按每初级归一；`verify.py` 仅用 `0.001` 将 `/GeV` 换算为 `/MeV`，乘精确能宽与立体角求和得 `0.0028799872629304804`，与 `fluka_unit_yields[0]` 一致，并输出 28 格明细表）。
- **Geant4 第 0 批原件**：[穿越记录 exits.csv](formal-g000-f00/geant4/exits.csv)、[初级记录 primaries.csv](formal-g000-f00/geant4/primaries.csv)（`10,000` 初级，窗口中 `49` 次中子出靶穿越）、[identity.tsv](formal-g000-f00/geant4/identity.tsv)、[incl-physics.tsv](formal-g000-f00/geant4/incl-physics.tsv)、[cuts-energy.tsv](formal-g000-f00/geant4/cuts-energy.tsv)、[rng_seed.txt](formal-g000-f00/geant4/rng_seed.txt) 与省略第 20 行集群路径后的[m1-readback 摘录](formal-g000-f00/geant4/m1-readback-excerpt.json)（其余共享输入文件与本目录 `geometry.gdml`/`beam.mac`/`physics.txt` 逐字节相同）。

单位换算（S07）与计分核查／归约的实现见[固定宿主](../settings-source/registered-main.cc)、[物理核查](../settings-source/base_kernel.py)、[计分核查](../settings-source/scoring.py)、[计分归约](../settings-source/score_reduction.py)、[USRYIELD 解析器](../settings-source/usryield.py)与[题面契约](../contracts/problem-fe.json)。开发期铁靶样本 Fe-M02 从 `20.534 mb` 截面型计分改为每初级产额型计分的回查与续跑记录见[早期计分案例](../records/early-Fe-scoring.json)。
