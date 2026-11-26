# Reading and writing a UMR mechanism

[Repository](../README.md) · [中文](language.zh-CN.md) · [Version 1.0 JSON Schema](../src/mechanism_workbench/resources/umr-v1.0.schema.json)

## Start with a Lewis structure

Open [substitution.umr](../examples/paper-1.0/substitution.umr). Its first line, `UMR 1.0 INITIAL s0;`, selects the language version and starting state. `STATE s0` declares the total charge and identifies `v0` as the primary view. `VIEW v0` lists an atom-mapped Lewis structure belonging to that state.

| Atom in `v0` | Attached H | Nonbonding electrons | Bond to another mapped atom |
|---|---:|---:|---|
| N1 of ammonia | 3 | 2 | — |
| C2 of methyl bromide | 3 | 0 | C2–Br3 |
| Br3 | 0 | 6 | Br3–C2 |

The map numbers persist through the event. Formal charge and nonbonding electron count are separate fields on each atom. Here all three atoms start neutral.

## Commit the arrows together

`REACT e1` binds N1, C2 and Br3 to the nucleophile, electrophile and leaving group of `SN2`. The following two lines are the action block in the complete example.

```text
PAIR(1,2);
SPLIT(2,3->3);
```

`PAIR` takes a lone pair from N1 to form N1–C2. `SPLIT` removes C2–Br3 and returns that pair to Br3. The checker considers the two arrows as a single event. Its successor `v1` has methylammonium (+1) and bromide (−1), with the same atom inventory and total charge as `v0`.

```sh
python umr.py examples/paper-1.0/substitution.umr --receipt sn2.json
python umr.py examples/paper-1.0/substitution-omitted-arrow.umr
```

The second input gives `Actions do not match SN2 combination rule`, returns exit code 1 and commits zero events. A model that receives this error can revise the event at the original intermediate. In a saved receipt, `document` contains the computed views and states, `event_reports` contains the checker results for each executed event and `errors` lists any execution errors. `full_electron_path_verified` reports whether the submitted path passed the strict electron-path checks. A declared net transformation can document a chemically observed reactant–product relation at a different level of detail.

## Read three other types of event

| Native input | How the arrows are represented | What to inspect in the product |
|---|---|---|
| [Radical substitution at sulfur](../examples/paper-1.0/radical-substitution.umr) | `FISHHOOK(RAD(1),NEW(1,3))` supplies one electron of the new C1–S3 bond. The two single electrons of S4–S3 have slots 0 and 1. One goes to the new bond and the other to S4. | The radical has moved from C1 to S4. |
| [Thermal Diels–Alder cycloaddition](../examples/paper-1.0/cycloaddition.umr) | Three `SHIFT(a,b,c)` actions each move a bonding pair from a–b to b–c. All three act in one six-electron event. | The product has new bonds 1–5 and 4–6 and a double bond between 2 and 3. |
| [Oxidative addition to Pd(0)](../examples/paper-1.0/oxidative-addition.umr) | `PAIR(20,1)` and `SHIFT(1,8,20)` record the formation of Pd–C and Pd–I as the aryl C–I bond breaks. | The state records Pd(II) and the metal–ligand connections. |

A change of Lewis representation has its own role. In the [enolate example](../examples/paper-1.0/enolate-branches.umr), the carbon- and oxygen-centred contributors are two views of the same state. A chemical event can start from either view and create a new state. For a site with an unusual electronic configuration, the [singlet silylene example](../examples/paper-1.0/silylene.umr) declares the Si centre's two bonds, one lone pair and zero unpaired electrons with `VALENCE`. Its `INTERMEDIATE` declaration tests both new Si–C bonds and the remaining C=C bond after capture of acetylene.

## Follow a named reaction through four events

In the [Larock program](../examples/paper-1.0/larock.umr), ortho-iodoaniline and 2-butyne form 2,3-dimethylindole. Starting from Pd(0), follow the aryl C–I oxidative addition, alkyne insertion into the aryl–Pd bond, acetate-assisted Pd–N bond formation, and C–N reductive elimination. The last event forms C12–N7 and regenerates Pd(0). The particular program combines alkyne coordination and insertion into one event.

```sh
python umr.py examples/paper-1.0/larock.umr --receipt larock.json --draw larock-diagram
python named_reactions.py --replay STR-260
```

The first command writes the structures and event reports to `larock.json`. The drawing directory contains `mechanism_chain.svg` and a separate SVG for each event. The second command checks the worked example attached to p. 260 in [Kürti and Czakó](../named-reactions/Kurti-Czako.md). For another family, try `python named_reactions.py --replay toc-002`. This replays the two declared events in the Aldol entry of [Li's catalogue](../named-reactions/Li.md). The textbook name identifies a reaction family; the receipt describes exactly what this particular substrate and program establish.
