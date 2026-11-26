# UMR 1.0

[English](README.md) · [语言入门](docs/language.zh-CN.md) · [人名反应库](#人名反应库)

**UMR：用于有机反应机理表示与结构验证的 Lean 式语言**  
邵思睿、谭海川

UMR 用带原子映射的 Lewis 结构和基元事件表示有机反应机理。化学家或语言模型提出电子如何移动，检查器核对电子来源、成断键、局部电子态、原子库存和总电荷。执行器据此生成下一中间体。同一组电子动作也给出了机理图中每根弯箭的起点和终点。

![从左向右依次为 UMR 指令、电子移动和后态，四行分别为极性取代、硫上的自由基取代、热环加成和氧化加成](docs/figures/elementary-events.zh.svg)

## 运行一个基元事件

克隆仓库后，在 Python 虚拟环境中从仓库根目录运行以下命令。若在无图形界面的服务器上绘图，可设置环境变量 `QT_QPA_PLATFORM=offscreen`。

```sh
python -m pip install -r requirements.txt
python umr.py examples/paper-1.0/substitution.umr
```

在这步 S<sub>N</sub>2 反应中，氨的 N1 用一对孤对电子与甲基碳 C2 成键，C2–Br3 键上的电子对归 Br3。两个动作在同一个事件内提交。

```text
PAIR(1,2);
SPLIT(2,3->3);
```

程序返回 `{"status": "passed", "events": 1, "strict": true, "errors": []}`。运行[删去第二根弯箭的输入](examples/paper-1.0/substitution-omitted-arrow.umr)，可以查看 S<sub>N</sub>2 动作组合检查的实际报错。要观察基元事件如何接成催化路径，可以沿 [Larock 吲哚合成](examples/paper-1.0/larock.umr)追踪钯的四次变化。

```sh
python umr.py examples/paper-1.0/larock.umr --receipt larock.json --draw larock-diagram
```

打开 `larock-diagram/mechanism_chain.svg` 即可查看连续机理图。`larock.json` 保存各步中间体和事件的执行结果。[语言入门](docs/language.zh-CN.md)用实际输入说明如何读回执、如何写自由基及周环事件，以及如何为单重态硅宾声明局部规则。

## 人名反应库

目录条目给出反应名称和教材页码。条目所链接的 UMR 程序进一步给出具体底物、原子映射和基元事件。

| 教材目录 | 收录范围 | 可回放内容 |
|---|---|---|
| [Jie Jack Li，《Name Reactions》第四版](named-reactions/Li.md) | 277 个主条目，另有 1 个单列的子条目 | 191 条可回放其声明的事件序列，80 条含部分路径或规则例外，1 条记录实测有机净转化，6 条仅列出处。 |
| [László Kürti、Barbara Czakó，《Strategic Applications of Named Reactions in Organic Synthesis》](named-reactions/Kurti-Czako.md) | 全部 250 章 | 48 个有具体底物的实例，共回放 210 个事件。 |

```sh
python named_reactions.py --find Larock
python named_reactions.py --show STR-260
python named_reactions.py --replay STR-260
```

Larock 条目对应第 260 页，回放结果为 `passed  events=4  strict=True`。程序核对输入文件哈希，并将本次状态、事件数和严格验证结果与条目记录逐项对照。`--check-all` 对全部已挂接程序进行同样的核对，包括保留部分路径和声明例外的记录。

## 与论文的对应关系

在 FukuyamaBench 已知终产物的机理推导任务中，模型获得起始物、终产物和有序条件，推导中间步骤。两组使用同一 GPT-6.2 Astra，四个示范的化学内容相同，每题各作答一次。**未使用 UMR 时命中全部官方有序检查点的有 62/319 题（19.4%），使用 UMR 与执行器时为 83/319 题（26.0%）。** 实验使用 UMR 0.9，本版语言示例使用 1.0。回放通过说明各步在声明的规则下结构自洽。哪条路径在给定条件下占优，仍须依据具体化学证据判断。

## 许可与引用

代码采用 MIT 许可。作者整理的人名反应目录、数据和原生程序采用 CC BY 4.0 许可，包括在 `src/` 和 `examples/` 中分发的副本。第三方组件沿用各自许可。文件范围和署名要求见[许可说明](NOTICE.md)，引用信息见 [CITATION.cff](CITATION.cff)。
