# SignatureAtlas

> **Release v4.** Versioning is by git tag (`v2`, `v3`, `v4`, …) — data files are unprefixed
> (`wg.bed`, `hm450.bed`, `defs.tsv`); check out a tag or a GitHub release for an earlier
> atlas. v4 (1,075 WGBS regions / 31,122 CpGs · 139 HM450 signatures / 11,760 probes)
> extends v3 from 88 to 100 cell-type/lineage contrasts; v3 introduced direct-clustering
> region selection, best-separation dedup, and the probe-span-seeded HM450 loci.

Whole-genome (WGBS) and HM450 DNA-methylation **marker signatures** for human cell
types and lineages, derived from an **87-cell-type** deep reference methylome. Each
signature is one **Most-Recurrent Methylation Pattern (MRMP)** — a recurrent binary
methylated/unmethylated partition of the panel — curated into compact genomic loci.
A signature is named by its **`signature_id`** = `<annotation><dir>`, where `dir` is
`-` (target **hypo**methylated) or `+` (target **hyper**methylated), e.g. `Squamous-`,
`CD8_T_Cell-`, `Oligodendrocyte+`.

## Workflow

```
   raw WGBS .cg stores (hg38, 29.4M CpGs)
   2023_Loyfer   2025_Zhou_MajorPseudo   2018_Zhou(BLUEPRINT)   2023_TianMajorType
        \________ per-cell-type / per-group musum pseudobulk ________/
                                    |
   [sa_00] --> sub.cg   87-cell reference panel        (yame = ~/repo/YAME/yame)
                                    |
 == Step A: MRMP layer =================================================
   [sa_01] keep CpG iff: all 87 covered & delta_mean>=0.6 & not chrY
           --> frequency-rank --> top-10,000 MRMPs (P1..P10000)
                                    |   each MRMP = one unique in_group/out_group partition
   [sa_02] + defs.tsv  (HAND-CURATED: in_group cells --> annotation name)
           exact-match each MRMP minority set --> signature_id = <annotation><+/->
           >> CONSISTENCY GATE: re-derived P-numbers must equal defs.tsv <<
                                    |
                  100 contrasts --> 169 signed signatures
                                    |
 == Step B: projection + curation ======================================
        +---------------------------+-------------------------------------+
        | HM450 (sparse array)      | whole-genome (WGBS, curated)
        v                           v
   [sa_03]  MRMP -> HM450 probes    [sa_04]  per-signature MRMP CpGs:
   [sa_05_hm450.sh] real beta         direct NR-distance clustering
   [sa_05_hm450.R]  keep a            (GAP = 100 CpGs, keep >=5-CpG runs)
       signature iff it has             + top-5 HM450-probe-dense loci (probe-span seeded)
       >= 5 probes                    --> base-curate on real beta (rowsub sub_small.cg)
        |                            [sa_refine_select]  worst-case outgroup separation
        |                            [sa_05_auto] region-gene (+/-10 kb) + sa_dedup_wg.R:
        |                              - MERGE overlapping same-signature regions
        |                              - BEST-SEPARATION: each CpG kept only in the
        |                                signature with the largest in-vs-out margin
        |                                (=> no duplicate coordinates)
        |                              - NAME each region by its first selected CpG
        |                            [sa_render_panels]  1 panel / region (merged windows)
        v                           v
   HM450/<signature_id>.ansi        WG/<signature_id>.<gene>_<chr>_<firstCpG>.ansi
   hm450.bed                     wg.bed
   139 signatures, 11,760 probes    164 signatures, 1,075 regions, 31,122 CpGs (0 dup)
```

## Layout

```
WG/<signature_id>.<gene>_<chr>_<firstCpG>.ansi   # per-region whole-genome panel
HM450/<signature_id>.ansi                        # per-signature HM450 panel
wg.bed                                         # all curated WGBS CpGs (1 row/CpG)
hm450.bed                                      # all selected HM450 probes (1 row/probe)
defs.tsv                                       # MRMP annotation (only hand-curated input)
```

- **`wg.bed`** — every curated whole-genome CpG, 5 columns: `chr`, `beg0`, `end1`
  (0-based begin / 1-based end; a CpG is the **2 bp** dinucleotide `(C-1, C+1)`, so
  `end1 - beg0 == 2`), `chr:min_max` (the region's full CpG span `firstCpG_lastCpG`),
  `region_name` (`<signature_id>.<gene>_<chr>_<firstCpG>`). **1,075** gene-named regions,
  **31,122** CpGs, **no duplicate coordinates**.
- **`hm450.bed`** — every selected HM450 probe, 5 columns: `chr`, `beg0`, `end1`
  (2 bp CpG), `probe_id` (cg number), `signature_id`. **11,760** probes across **139**
  signatures.
- **`WG/*.ansi` / `HM450/*.ansi`** — colored real-β inspection panels (**panels only —
  no per-signature beds**): stacked HM450-probe track (cg_id-labelled) + the selection
  (SEL) track + the 87-cell β heatmap, target rows **underlined**, cropped to the
  selection ± 50 flanking CpGs. WG panels are **1:1** with `wg.bed` regions over the
  merged β window; a CpG reassigned away by best-separation still appears in the donor
  panel's heatmap (context) but not in its SEL track.

## Coordinate & naming conventions

- **2 bp CpG, 0-based beg / 1-based end** in both BEDs: for a CpG whose C is at 1-based
  `pos`, `beg0 = pos - 1`, `end1 = pos + 1` (`end1 - beg0 == 2`).
- **Region identity = its selected CpGs.** A WGBS region is named by its **first
  (lowest-coordinate) selected CpG** plus the dominant gene within ±10 kb; the
  `chr:min_max` span is `firstCpG_lastCpG`. The 10 kb HM450 bucketing only *ranks*
  array-dense loci to seed extraction — it never names a region.
- **No duplicate coordinates.** Overlapping regions of the **same** signature are
  merged into one; a CpG selected by **multiple** signatures is kept only in its
  **best-separation** signature (largest worst-case in-vs-out β margin) and dropped
  from the others.

## defs.tsv

The **MRMP annotation** — the only hand-curated input to the atlas. Each row defines a
**contrast** (an `in_group`/`out_group` partition of the panel) and gives it a
hand-chosen `annotation` name; **everything else is automatic** (matching MRMPs, HM450
probe / whole-genome region selection, curation, refine). One row per contrast; columns:
`annotation`, `n_in`, `mrmp_neg`, `mrmp_pos`, `out_group`, `in_group`.

- `annotation`, `in_group` — *(curated)* the contrast name and the comma-list of panel
  cell types on its `in_group` side. Cell-type **equivalence is implicit** in the list
  (e.g. `NK_cell,Hema_NK,BP_NKCell`), and a contrast may be a single cell, a lineage, or
  any hand-chosen set (vascular, microglia, inhibitory-neuron, oligodendrocyte-lineage…).
- `mrmp_neg` / `mrmp_pos` — *(automatic)* the P-numbers of the MRMPs realizing this
  contrast as a `-` (in_group hypo) or `+` (in_group hyper) **signature**, i.e. the MRMP
  whose minority set **equals** `in_group`. A contrast matching no MRMP yields no
  signature and is dropped.
- `n_in`, `out_group` — *(automatic)* the in_group size and the out-group (`COMPLEMENT`
  = every other cell).

Currently **100 contrasts → 169 signatures** (each realized `+`/`-` direction is a signed
`signature_id`); **164** are realized as curated WGBS regions and **139** carry ≥5 HM450
probes. Only the curated `defs.tsv` + the panel (`sub.cg`) drive everything
downstream; edit a contrast's `annotation`/`in_group`, then re-run from `sa_02`. The
`sa_02` consistency gate re-derives the P-numbers from `in_group` and **fails the build**
if they disagree with the `mrmp_neg`/`mrmp_pos` already in `defs.tsv` — so the
hand-curated P-numbers stay locked to the panel.

## Methods (summary)

**Reference panel (87 cell types).** 31 Loyfer per-cell-type pseudobulks
(`yame rowop -o musum`) + 34 Zhou "Major" pseudobulks + **8 BLUEPRINT** disease-free
immune pseudobulks (Monocyte, Macrophage, Dendritic, Neutrophil, CD4/CD8 T, NK, B) +
**14 Tian** brain groups (excitatory + Pvalb/Sst/Vip/Lamp5/MSN interneurons + ASC/ODC/
OPC/MGC/EC/PC/VLMC glia & vascular). The BLUEPRINT + Tian additions give immune and
brain granularity beyond the Loyfer/Zhou base.

**MRMP discovery (`sa_01`).** At each of 29,401,795 CpGs, keep sites with all 87 cell
types covered, high-vs-low group-mean Δ ≥ 0.6 (`delta_mean`, a *derivation* filter —
every MRMP CpG inherits it; emitted by `yame rowop -o stat`), and not on chrY; reduce to
an 87-bit methylated/unmethylated string; rank by genome-wide frequency; keep the top
10,000. Each P-number is just a unique `in_group`/`out_group` partition; names are
applied only in `sa_02`.

**Annotation (`sa_02`).** An MRMP already *defines* its contrast; annotation gives it a
biological **name** by matching the MRMP minority set **exactly** against the `in_group`
lists in `defs.tsv`, assigning `signature_id = <annotation><direction>`. HM450 and
whole-genome share the *same* signatures, differing only in projection (sparse array vs
curated genome-wide).

**Whole-genome curation (`sa_04` → `sa_05_auto`).** For each signature, its own MRMP
CpGs are clustered **directly** in genomic order — consecutive same-signature CpGs within
**GAP = 100 CpG-indices** join a run; runs of ≥5 CpGs are candidate DMR seeds (top-5 by
size) — plus the top-5 HM450-probe-dense 10 kb loci (so the WGBS view also covers
array-anchored regions). Each seed is base-curated on **real β** (`yame rowsub` of a
pooled `sub_small.cg` → `unpack -a -f 5`; hypo β<0.4 in the in_group / hyper β≥0.6) with a
±100-CpG flank, then `sa_refine_select` keeps only CpGs whose worst-case smoothed
out-group separation holds. `sa_dedup_wg.R` then merges overlapping same-signature
regions, resolves cross-signature CpGs by best separation, and names each region by its
first selected CpG. (Win10k window *enrichment* — used in earlier versions to pick
regions — has been removed in favour of this direct clustering.)

**HM450 projection (`sa_03` → `sa_05_hm450`).** MRMP CpGs that carry an HM450 probe are
taken per signature; a signature is kept iff it retains ≥5 probes. (Every MRMP CpG already
separates in_group from out_group by `delta_mean` ≥ 0.6, so no extra per-probe
effect-size filter is needed.)

## Pipeline scripts (in order)

```
Step A:  sa_01_mrmp.sh -> sa_01_def.R -> sa_02_classify.R   (+ consistency gate)
Step B:  sa_03_hm450_win.sh -> sa_05_hm450.sh -> sa_05_hm450.R
         sa_04_enrich_curate.sh -> sa_refine_select.R
         sa_05_auto.sh (region_gene + sa_dedup_wg.R) -> sa_render_panels.R
```

## TODO / Known issues

- **Curate the unannotated MRMPs in `mrmp_annotated.tsv`.** Only 169 of the 10,000 MRMPs
  match a curated contrast in `defs.tsv`; the rest are unnamed partitions (`signature_id`
  = `NA`). Adding contrasts to `defs.tsv` for the recurrent unannotated patterns would
  extend cell-type/lineage coverage.
- **`skeletal_muscle` (Loyfer) appears impure** (sorting / whole-tissue): no standalone
  MRMP, only co-segregates with Zhou `Mus_Skl` ± scattered epithelial contaminants. The
  exact-partition MRMP scheme splits the skeletal signal (`{Mus_Skl}` only vs
  `{Mus_Skl, skeletal_muscle}`), so the signature loses sites. **Fix (deferred):** drop
  `skeletal_muscle` and re-derive so `Mus_Skl` consolidates.
- **`cortex_neuron` (Loyfer) is impure** and biases the broad `Neuronal` contrast:
  no standalone MRMP, excitatory-dominated, drags in `Fibro_Brst` contamination.
  **Fix (deferred):** exclude it from the `Neuronal` signature (re-derive without it).
- **`OPC+` has no WGBS region** (rank-8228 MRMP, scattered CpGs form no ≥5-CpG cluster
  and it carries no HM450 probe). Accepted — no robust local DMR exists for it.
