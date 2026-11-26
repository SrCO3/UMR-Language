# UMR 1.0

[简体中文](README.zh-CN.md) · [Language walkthrough](docs/language.md) · [Named reaction library](#named-reaction-library)

**UMR: A Lean-Style Language for Representation and Structural Verification of Organic Reaction Mechanisms**  
Shao Sirui and Tan Haichuan

UMR represents a mechanism as a sequence of elementary events on atom-mapped Lewis structures. A chemist or language model proposes how electrons move. The checker tests the electron sources, bond changes, local electronic configurations, atom inventory and total charge. The executor constructs the next intermediate. The same actions specify the arrows in a mechanism diagram.

![UMR instructions, electron flow and resulting structures for polar substitution, radical substitution, cycloaddition and oxidative addition](docs/figures/elementary-events.en.svg)

## Run an elementary event

Clone this repository and run the commands in its root using Python in a virtual environment. The drawing command below also works on a headless machine with `QT_QPA_PLATFORM=offscreen`.

```sh
python -m pip install -r requirements.txt
python umr.py examples/paper-1.0/substitution.umr
```

In this S<sub>N</sub>2 reaction, ammonia donates its lone pair on N1 to methyl carbon C2, and the C2–Br3 bonding pair moves to Br3. The two actions commit together.

```text
PAIR(1,2);
SPLIT(2,3->3);
```

The program returns `{"status": "passed", "events": 1, "strict": true, "errors": []}`. Try the [input with the second arrow omitted](examples/paper-1.0/substitution-omitted-arrow.umr) to see the rule that rejects the incomplete event. For a four-event example, follow Pd through the [Larock indole synthesis](examples/paper-1.0/larock.umr).

```sh
python umr.py examples/paper-1.0/larock.umr --receipt larock.json --draw larock-diagram
```

Open `larock-diagram/mechanism_chain.svg` for the complete path. The JSON receipt records the intermediate structures and the outcome of each event. The [walkthrough](docs/language.md) explains how to read a receipt, what changes in a radical or pericyclic event, and how to declare a local rule for a singlet silylene.

## Named reaction library

Each catalogue entry identifies a reaction and its printed page. A linked UMR program specifies the substrate, atom mapping and events for a particular example.

| Catalogue | Entries | Executable material |
|---|---|---|
| [Jie Jack Li, *Name Reactions*, 4th ed.](named-reactions/Li.md) | 277 main entries and one separately indexed subentry | 191 records replay their declared sequence, 80 contain a partial path or an exception, one describes the measured organic net transformation, and six provide a source reference. |
| [László Kürti and Barbara Czakó, *Strategic Applications of Named Reactions in Organic Synthesis*](named-reactions/Kurti-Czako.md) | 250 printed chapters | 48 substrate-specific examples with 210 replayed events. |

```sh
python named_reactions.py --find Larock
python named_reactions.py --show STR-260
python named_reactions.py --replay STR-260
```

The Larock entry is on p. 260. Replay gives `passed  events=4  strict=True`. The tool verifies the input digest and compares the status, event count and strictness with the catalogue record. `--check-all` applies this comparison to every linked program, including records with partial paths and declared exceptions.

## Relation to the paper

In FukuyamaBench *Pathway-only Prediction*, the model receives the starting materials, known product and ordered conditions, then proposes the intervening steps. With the same four chemical demonstrations and one attempt per question, GPT-6.2 Astra matched all official ordered checkpoints on **62/319 questions (19.4%) without UMR** and **83/319 (26.0%) using UMR and the executor**. The experiment used UMR 0.9. The language examples in this release use 1.0. A successful replay establishes structural consistency under the declared rules. Selecting the pathway favored under the experimental conditions also requires chemical evidence.

## License and attribution

The code is licensed under MIT. The authors' named-reaction catalogue, data and programs are licensed under CC BY 4.0, including copies distributed under `src/` and `examples/`. Third-party components retain their own licenses. See the [notice and file-level scope](NOTICE.md) and [citation](CITATION.cff).
