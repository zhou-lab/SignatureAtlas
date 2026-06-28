# SignatureAtlas

Whole-genome (WGBS) and HM450 DNA-methylation **marker signatures** for human cell
types and lineages, derived from an **87-cell-type** deep reference methylome. Each
signature is one **Most-Recurrent Methylation Pattern (MRMP)** — a recurrent
binary methylated/unmethylated partition of the panel — curated into compact genomic
loci. A signature is named by its **`signature_id`** = `<annotation><dir>`, where
`dir` is `-` (target **hypo**methylated) or `+` (target **hyper**methylated), e.g.
`Squamous-`, `CD8_T_Cell-`, `Oligodendrocyte+`.

## Workflow

```
   raw WGBS .cg stores (hg38, 29.4M CpGs)
   2023_Loyfer  2025_Zhou_MajorPseudo  2018_Zhou(BLUEPRINT)  2023_TianMajorType
        \________ per-cell-type / per-group musum pseudobulk ________/
                                   |
   [sa_00] -> sub.cg   87-cell reference panel
                                   |
   [sa_01] keep CpG iff: all 87 covered & delta_mean>=0.6 & not chrY
           -> frequency-rank -> top-10,000 MRMPs
                                   |   each MRMP (P#) = one unique in_group/out_group partition
   [sa_02] + v2_defs.tsv  (HAND-CURATED: in_group cells -> signature name)
           exact-match each MRMP in_group -> signature_id = <annotation><+/->
                                   |
        +--------------------------+----------------------------------+
        | HM450 (sparse array)     | whole-genome (WGBS, curated)
        v                          v
   [sa_03 -> sa_05_hm450]     [sa_04 enrich+curate -> refine -> sa_05_auto]
     per-probe delta-beta       Win10k enrich -> top-5 windows -> decile base-curation
        v                       -> worst-case specificity refine -> per-region gene
   HM450/<sig>.ansi              v
   v2_hm450.bed                  WG/<sig>.<gene>_<chr>_<beg>.ansi
                                 v2_wg.bed
```

## Layout

```
HM450/<signature_id>.ansi                 # per-signature HM450 inspection panel
WG/<signature_id>.<gene>_<chr>_<beg>.ansi # per-region whole-genome inspection panel
v2_wg.bed                                 # all curated whole-genome CpGs (1 row/CpG)
v2_hm450.bed                              # all selected HM450 probes (1 row/probe)
v2_defs.tsv                               # the MRMP annotation (only hand-curated input)
ARCHIVE/                                  # frozen previous atlases (flattened beds)
  v0_hm450.bed                            #   v0 (2026-02-03) HM450 signatures
  v1_hm450.bed                            #   v1 (65-cell) HM450
  v1_wg.bed                               #   v1 (65-cell) whole-genome
```

- **`v2_wg.bed`** — every curated whole-genome CpG, 5 columns:
  `chr`, `beg0`, `end1` (0-based begin / 1-based end; single CpG = `(C-1, C)`),
  `chr:min_max` (the full span of that signature's region), `signature` (the region
  name `<signature_id>.<gene>_<chr>_<beg>`). 582 gene-annotated regions.
- **`v2_hm450.bed`** — every selected HM450 probe, 5 columns: `chr`, `beg0`,
  `end1`, `probe_id` (cg number), `signature_id`. 11,105 probes across 112 signatures.
- **`HM450/*.ansi` / `WG/*.ansi`** — colored decile/β inspection panels: stacked
  HM450-probe track (cg_id labelled) + the selection (SEL) track + the 87-cell
  heatmap, target rows **underlined**, cropped to the selection ± 50 flanking CpGs.

## v2_defs.tsv

The **MRMP annotation** — the only hand-curated input to the atlas. Each row defines a
**contrast** (an `in_group`/`out_group` partition of the panel) and gives it a
hand-chosen `annotation` name; **everything else is automatic** (matching MRMPs, HM450
probe / whole-genome region selection, curation, refine). One row per contrast; columns:
`annotation`, `n_in`, `mrmp_neg`, `mrmp_pos`, `out_group`, `in_group`.

- `annotation`, `in_group` — *(curated)* the contrast name and the comma-list of panel
  cell types on its `in_group` side. Cell-type **equivalence is implicit** in the list
  (e.g. `NK_cell,Hema_NK,BP_NKCell`), and a contrast may be a single cell, a lineage, or
  any hand-chosen set (vascular, microglia, inhibitory-neuron, oligodendrocyte-lineage, …).
- `mrmp_neg` / `mrmp_pos` — *(automatic)* the P-numbers of the MRMPs realizing this
  contrast as a `-` (in_group hypo) or `+` (in_group hyper) **signature**, i.e. the MRMP
  whose minority **equals** `in_group`. A contrast matching no MRMP yields no signature
  and is dropped.
- `n_in`, `out_group` — *(automatic)* the in_group size and the out-group (`COMPLEMENT`
  = every other cell).

Currently **83 contrasts → 144 signatures** (each realized `+`/`-` direction is a signed
`signature_id`). Only the curated `v2_defs.tsv` + the panel (`sub.cg`) drive
everything downstream; edit a contrast's `annotation`/`in_group`, then re-run from `sa_02`.

## Methods (summary)

**Reference panel (87 cell types).** 31 Loyfer per-cell-type pseudobulks
(`yame rowop -o musum`) + 34 Zhou "Major" pseudobulks + **8 BLUEPRINT** disease-free
immune pseudobulks (Monocyte, Macrophage, Dendritic, Neutrophil, CD4/CD8 T, NK, B) +
**14 Tian** brain groups (excitatory + Pvalb/Sst/Vip/Lamp5/MSN interneurons + ASC/ODC/
OPC/MGC/EC/PC/VLMC glia & vascular). The BLUEPRINT + Tian additions give immune and
brain granularity beyond the Loyfer/Zhou base.

**MRMP discovery.** At each of 29,401,795 CpGs, keep sites with all 87 cell types
covered, high-vs-low group-mean Δ ≥ 0.6 (`delta_mean`, a *derivation* filter — every
MRMP CpG inherits it), and not on chrY; reduce to an 87-bit methylated/unmethylated
string; rank by genome-wide frequency; keep the top 10,000. The MRMP set carries **no**
label — each P-number is just a unique `in_group`/`out_group` partition; names are
applied only in `sa_02`.

**Annotation.** An MRMP already *defines* its contrast (the `in_group`/`out_group`
partition); annotation just gives that partition a biological **name**. Each MRMP's
minority set is matched **exactly** against the `in_group` lists in `v2_defs.tsv`;
a match assigns `signature_id = <annotation><direction>`. One MRMP ⇒ one signed
signature; HM450 and whole-genome share the *same* signatures, differing only in
projection (sparse array vs curated genome-wide).

## TODO / Known issues

- **`skeletal_muscle` (Loyfer) appears impure** (sorting / whole-tissue): no standalone
  MRMP, only co-segregates with Zhou `Mus_Skl` ± scattered epithelial contaminants. The
  exact-partition MRMP scheme splits the skeletal signal (`{Mus_Skl}` only vs
  `{Mus_Skl, skeletal_muscle}`), so the signature loses ~3–4k sites. **Fix (deferred):**
  drop `skeletal_muscle` and re-derive so `Mus_Skl` consolidates. Low priority — the
  both-hypo signature is adequate.
- **`cortex_neuron` (Loyfer) is impure** and biases the broad `Neuronal` contrast:
  no standalone MRMP, excitatory-dominated, drags in `Fibro_Brst` contamination.
  **Fix (deferred):** exclude it from the `Neuronal` signature (re-derive without it).
