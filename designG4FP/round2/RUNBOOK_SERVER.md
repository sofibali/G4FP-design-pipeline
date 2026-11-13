# Round-2 runbook — lab server (Sofia)

Everything in this runbook runs where the AlphaFold3 predictions live. Nothing
here needs a GPU: it reads the CIFs you already have and does geometry and
bookkeeping.

**Why the division of labour changed.** All 1,000 round-1 AF3 predictions are
on the lab server. The Drive export has only 36 designs (3 per template, 9 for
`r20_l1_cro_mod0`) and this repo has never contained a single `.cif` —
predictions are excluded by `designG4FP/.gitignore`. So the rescoring and
selection have to run where the structures are; only the small outputs come
back to git.

**What this replaces.** The round-1 chromophore pLDDT is wrong (see README §1.1:
it averages residues ~201-204 instead of the chromophore, overestimating by
12.8 points on average). Re-running the steps below regenerates correct values
for all 1,000 designs and selects round-2 inputs from them.

---

## 0. Prerequisites

```bash
conda activate unified_mpnn          # or any env with numpy + pandas
cd /home/sbali/LigandMPNN/designG4FP
git fetch origin && git checkout round2-redesign
cd round2
python -c "import numpy, pandas; print('ok')"
```

No other dependencies. The scripts parse mmCIF directly rather than via
BioPython, so there is nothing to install.

Set these two to match your filesystem:

```bash
export HOLO=/path/to/output_G4FP_*/03_alphafold3_predictions_bound
export APO=/path/to/output_G4FP_*/03_alphafold3_predictions_apo
```

If the predictions sit under one parent, just point at the parent — discovery
walks the tree and reads the template name from the enclosing `output_G4FP_*`
directory. Both the lab layout
(`output_G4FP_<tpl>/..._predictions_<state>/design_<NNNN>/seed-N_sample-M/*model*.cif`)
and a flat `<tpl>/<prefix>_d<NNNN>_s<k>.cif` export are recognised.

---

## 1. Rescore every design with corrected chromophore metrics

```bash
export TPL_APO=/path/to/template_af3_outputs      # apo predictions of the 10 TEMPLATES

python 09_rescore_and_select.py \
    --af3-root          "$HOLO" \
    --apo-root          "$APO" \
    --template-apo-root "$TPL_APO" \
    --layout            af3 \
    --round1-table      /path/to/Drive/07_all_designs_scored.csv \
    --validate          tested_variants.csv \
    --out               09_results_full \
    --top-n             5
```

**`--template-apo-root` is not optional.** `apo_plddt_vs_template` is what
drives the cascade's foldability gate, and it is **NaN for all 1,000 designs in
the committed `../07_results/plot_data_csvs/07_all_designs_scored.csv`** —
that is round-1 defect (3) in README §1.3, not something this pipeline has
fixed by itself. The script now computes the metric directly from structures:
each design's apo mean pLDDT minus its parent template's. That needs AF3 apo
predictions of the ten templates, which is what `$TPL_APO` points at.

Two fallbacks if you cannot produce template apo predictions: point
`--round1-table` at the **Drive** copy of `07_all_designs_scored.csv`
(1,020 rows, metric populated) rather than the committed one, or re-predict the
ten templates apo. If neither source supplies it the script now prints a
block-capital warning and the gate will pass everything — do not proceed past
that warning.

Expect roughly 25 structures per design (5 seeds x 5 samples) x 1,000 designs
x 2 states. Runtime is dominated by file I/O; budget an hour or two and run it
under `nohup`.

**Check before moving on**

- `Collected metrics for N designs across 10 templates` — N should be ~1000,
  not 36. If it says 36 you have pointed at the Drive export.
- `computed apo_plddt_vs_template from structures for N/N designs, against 10
  template references` — if you instead see `*** apo_plddt_vs_template is
  EMPTY for every design ***`, stop: stage B cannot run and the rel_apo term is
  inert. Fix the inputs before continuing.
- `Merged ... (N/N designs matched)` — if the match count is low, the
  `global_id` convention differs; `norm_id()` in the script strips the `G4FP_`
  prefix and zero-pads the index, so check what your directories are called.
- In `09_results_full/09_all_designs_rescored.csv`, `chromophore_start` should
  be 187, 195, 196 or 197 depending on the template — **not** uniformly 197.
  That column is how you know the motif finder located the real triad.
- Compare `holo_chromophore_plddt_fixed` against the old
  `holo_chromophore_plddt`: the corrected values should be systematically
  lower, by ~6-20 points.

**The one result worth looking at immediately.** The validation block prints
how the 21 experimentally tested variants rank. Round 1 ranked the two
fluorescent ones 15th and 19th of 21. If the corrected metrics still rank them
near the bottom, that is informative and expected — see README §2.

---

## 2. Filter with the iterative pipeline

```bash
python 12_iterate.py \
    --table     09_results_full/09_all_designs_rescored.csv \
    --iteration 1 \
    --validate  tested_variants.csv \
    --out       iter1
```

Four gate groups — APO (a window), FOLD, BIND, SWITCH — with percentile
thresholds that tighten at iteration 2 and 3. The per-criterion trace prints
as it runs and is written to `iter1/12_gate_trace.csv`.

On later rounds add `--previous iter1/12_survivors.csv`, which enforces that a
design passed **both before and after** the design/prediction cycle.

**Check before moving on**

- The APO floor should remove roughly 70% of the pool (730 of 1,010 in the
  demo run). If it removes nothing, `apo_plddt_vs_template` is NaN — go back
  to step 1.
- Every group should print real thresholds, not `metric unavailable --
  SKIPPED`. A skipped SWITCH group means no ligand-dependence is being tested
  at all, which defeats the purpose; it needs `--apo-root` in step 1.
- Expect of order 10-50 survivors at iteration 1 from a 1,000-design pool.

`iter1/12_next_round_inputs.csv` is the handoff to step 3.

### 2b. (superseded) `11_cascade.py`

Kept for reference only. It has no SWITCH group, so a construct equally folded
with and without G4 passes every gate — the pKZ168 phenotype. Use
`12_iterate.py`.

## 2c. Select round-2 inputs with the staged cascade (legacy)

```bash
python 11_cascade.py \
    --table    09_results_full/09_all_designs_rescored.csv \
    --validate tested_variants.csv \
    --out      11_results_full \
    --select   stratify \
    --rank-by  apo_plddt_vs_template
```

Gates are AND-ed and non-compensatory; the per-stage survivor counts print as
it runs. The chromophore gate (stage A) is off by default — see README §5.1 for
why. Every threshold is a flag if you want to explore:

```bash
--apo-vs-template 0.0     # stage B: apo at least as ordered as the template
--holo-plddt 70 --holo-ptm 0.5 --holo-iptm 0.3
--interface-pae 5         # only applies if the column is present
--per-band 2 --strat-edges 0,0.5,1,2,100
```

**Check before moving on**

- On the 36 designs I could test locally: 36 → 8 survivors, containing both
  fluorescent variants and none of the twelve tested dark ones.
- Stage B should remove most of the pool. 280 of 1,000 clear
  `apo_plddt_vs_template >= 0` **in the Drive copy** of the round-1 table, so a
  few hundred survivors is the right order of magnitude. If stage B removes
  **nothing**, the metric is all-NaN and the gate is silently passing
  everything — go back to step 1.
- If a stage prints `column absent -- criterion SKIPPED`, that criterion did
  not run. Interface PAE in particular is only present if your AF3 summary
  JSONs sit next to the model CIFs.
- If a stage drops everything, loosen that stage's flag rather than the others.

---

## 3. Graft the mature chromophore into the selected designs

```bash
python 10_graft_chromophore.py \
    --models    "$HOLO" \
    --templates ../inputs \
    --only      iter1/12_next_round_inputs.csv \
    --out       inputs_round2_full
```

AF3 does not model the mature chromophore — it predicts Thr-Tyr-Gly as an
extended tripeptide (Thr-Ca to Gly-Ca 4.66-5.35 A against 3.10 A in the
cyclized CRO). This transplants CRO from each design's parent template and
writes PDB, which LigandMPNN needs anyway since `run.py` reads PDB only.

**Check before moving on**

- `graft_fit_rmsd` should be well under 1 A for most designs. Anything above
  ~1.5 A means the barrel has moved relative to its template and the placement
  is not trustworthy — exclude those rather than designing on them.
- `new_clashes` at *good* fit is a real signal: the `r20_l2_cro_mod3` designs
  clash 11-15 times at fit 0.25-0.31 A, meaning that barrel genuinely cannot
  host a mature chromophore. A clash count at poor fit means nothing.

---

## 4. Commit the outputs back

The grafted PDBs are small and are the actual handoff:

```bash
git add round2/09_results_full/*.csv \
        round2/11_results_full/*.csv \
        round2/inputs_round2_full/
git commit -m "Round-2 rescoring and selection over the full prediction set"
git push origin round2-redesign
```

Do **not** commit the AF3 predictions themselves — `.gitignore` excludes them
deliberately and they are tens of GB.

---

## 5. Also needed: seven missing structures

Seven of the 21 experimentally tested variants have no AF3 structure in Drive
or anywhere in git history. They are listed by exact design id in
`MISSING_STRUCTURES.md`. If they still exist on the server, they fold into
step 1 automatically; if not, please re-predict them in the holo state
(protein + 22-nt Myc G4 + 3 K+, 5 seeds, round-1 settings).

---

## Known gotchas

| Symptom | Cause / fix |
|---|---|
| `Collected metrics for 36 designs` | pointed at the Drive export, not the server tree |
| template recorded as `03_alphafold3_predictions_bound` | fixed in commit on this branch — make sure you are on `round2-redesign` at or after `10b38c5` |
| `chromophore_start` uniformly 197 | the `VTTL(TYG)VQ` motif did not match; check the construct sequence, then use `--chromophore-fallback 197,199` only as a last resort |
| `rel_apo term is inert` / `apo_plddt_vs_template is EMPTY` | the metric is unavailable. Pass `--apo-root` **and** `--template-apo-root` to compute it from structures, or point `--round1-table` at the Drive copy. The committed `07_all_designs_scored.csv` has it all-NaN (README §1.3) |
| stage B removes zero designs | same cause: the gate cannot evaluate a NaN column and passes everything |
| all designs removed by a gate | loosen that gate's flag; the printout names which criterion removed what |
| `No structures found under ...` | wrong `--layout`; try `--layout af3` explicitly |

## IMPORTANT: --apo-root is not optional

`09_rescore_and_select.py` must be given `--apo-root` pointing at the matching
apo predictions. Without it `chromophore_plddt_diff_fixed` and
`chromophore_rmsd` are never computed, the whole SWITCH gate is silently
skipped, and `12_iterate.py` degenerates to APO + FOLD + ipTM. The script
prints `metrics unavailable (their criteria are skipped): ...` when this
happens - check that line before trusting a survivor list.

`--template-apo-root` is likewise required for `apo_plddt_vs_template`.
