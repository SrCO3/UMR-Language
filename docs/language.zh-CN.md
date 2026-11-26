# 阅读与编写 UMR 机理

[仓库首页](../README.zh-CN.md) · [English](language.md) · [1.0 版 JSON Schema](../src/mechanism_workbench/resources/umr-v1.0.schema.json)

## 从 Lewis 结构写起

打开 [substitution.umr](../examples/paper-1.0/substitution.umr)。第一行 `UMR 1.0 INITIAL s0;` 指定语言版号和初态。`STATE s0` 声明体系总电荷，选定 `v0` 为该状态的主要 Lewis 表示。`VIEW v0` 写出属于这个状态的原子和键。

| `v0` 中的原子 | 所连氢数 | 非键电子数 | 与其他映射原子之间的键 |
|---|---:|---:|---|
| 氨中的 N1 | 3 | 2 | 无 |
| 溴甲烷中的 C2 | 3 | 0 | C2–Br3 |
| Br3 | 0 | 6 | Br3–C2 |

同一原子在反应前后保留相同编号。每个原子的形式电荷与非键电子数分别记录。这里三个原子的初始形式电荷都是零。

## 将同一步的弯箭一起提交

`REACT e1` 指定事件类型 `SN2`，并把 N1、C2、Br3 分别绑定到亲核位点、亲电位点和离去基团。完整输入中的动作如下。

```text
PAIR(1,2);
SPLIT(2,3->3);
```

`PAIR` 将 N1 的一对孤对电子用于形成 N1–C2 键。`SPLIT` 断开 C2–Br3，键电子归 Br3。检查器将两根弯箭作为一个事件处理。后态 `v1` 包含甲铵离子（+1）与溴离子（−1），原子库存和总电荷均与初态一致。

```sh
python umr.py examples/paper-1.0/substitution.umr --receipt sn2.json
python umr.py examples/paper-1.0/substitution-omitted-arrow.umr
```

第二个文件删除了离去基团对应的弯箭。程序输出 `Actions do not match SN2 combination rule`，退出码为 1，已执行事件数为 0。语言模型收到这个错误后，可以从原中间体修改该事件。回执中的 `document` 保存执行得到的各状态和 Lewis 表示，`event_reports` 保存逐事件检查结果，`errors` 收录执行错误。`full_electron_path_verified` 表明所提交的电子路径是否整体通过严格检查。对于实验上观察到的底物与产物关系，还可以按净转化记录，并注明机理展开到了哪一层。

## 阅读另外三类基元事件

| 原生输入 | 电子移动怎样写 | 产物中核对什么 |
|---|---|---|
| [硫上的自由基取代](../examples/paper-1.0/radical-substitution.umr) | `FISHHOOK(RAD(1),NEW(1,3))` 把 C1 的未成对电子送到新 C1–S3 键。原 S4–S3 键上的两个电子分别编号为 0 和 1，一个进入新键，另一个留在 S4。 | 自由基中心由 C1 转移至 S4。 |
| [热 Diels–Alder 环加成](../examples/paper-1.0/cycloaddition.umr) | 每条 `SHIFT(a,b,c)` 将 a–b 键上的一对电子移到 b–c。三条动作组成同一六电子事件。 | 生成 1–5、4–6 两根新键，2 与 3 之间形成双键。 |
| [Pd(0) 氧化加成](../examples/paper-1.0/oxidative-addition.umr) | `PAIR(20,1)` 和 `SHIFT(1,8,20)` 共同描述芳基 C–I 键断裂、Pd–C 与 Pd–I 键形成。 | 后态给出 Pd(II) 和金属与配体的连接。 |

同一物种还可以有多个 Lewis 表示。[烯醇负离子实例](../examples/paper-1.0/enolate-branches.umr)把负电荷分别写在碳、氧上的两个共振式归在同一状态下。化学事件可以从其中任一表示出发，生成新的状态。遇到特殊电子态时，可以按位点声明规则。[单重态硅宾实例](../examples/paper-1.0/silylene.umr)用 `VALENCE` 指定 Si 连两根键、有一对孤对电子、未成对电子数为零。捕获乙炔后，`INTERMEDIATE` 核对两根新 Si–C 键和保留的 C=C 键。

## 沿四个事件读人名反应

在 [Larock 程序](../examples/paper-1.0/larock.umr)中，邻碘苯胺与 2-丁炔形成 2,3-二甲基吲哚。沿 Pd(0) 依次追踪芳基碘的氧化加成、炔烃向芳基–Pd 键的插入、乙酸根参与的 Pd–N 成键及 C–N 还原消除。最后一步形成 C12–N7 键，再生 Pd(0)。该实例将炔烃配位和迁移插入合写为一个事件。

```sh
python umr.py examples/paper-1.0/larock.umr --receipt larock.json --draw larock-diagram
python named_reactions.py --replay STR-260
```

第一条命令将各状态结构及逐事件回执写入 `larock.json`。绘图目录中的 `mechanism_chain.svg` 是连续机理图，另有每个事件独立的 SVG。第二条命令核对 [Kürti/Czakó 目录](../named-reactions/Kurti-Czako.md)第 260 页所挂接的具体实例。也可以运行 `python named_reactions.py --replay toc-002`，回放 [Li 目录](../named-reactions/Li.md)中 Aldol 条目声明的两个事件。教材名称指向反应家族，回执则把所用底物、执行了哪些步、哪些检查通过落实到具体实例。
