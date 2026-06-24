# SignatureAtlas

Whole-genome and HM450 DNA-methylation **marker signatures** for human cell types,
lineages, and within-lineage subgroups, derived from a **65-cell-type** deep WGBS
reference (31 Loyfer per-cell-type pseudobulks + 34 Zhou "Major" pseudobulks across
five lineages: Epith / Hema / Neural / Mesench / Adrenal). Markers are
Most-Recurrent Methylation Patterns (MRMPs) curated into compact genomic regions,
labelled `<cell_type>` (cell-type-specific), `pan<Lineage>` (pan-lineage), or
`sub<Lineage>_<sublineage>` (within-lineage subgroup), with direction `-` =
hypomethylated and `+` = hypermethylated in the target.

The full reproducible pipeline (control sheet → this atlas) lives in the lab
journal: `labjournal/zhouw3/2026/20260619_signature_atlas.org` + `tools/sa_0*`.

## Layout

```
auto/<cell_type|panLineage|subLineage_sublineage>/<label><dir>.<gene>_<chr>_<beg>.{bed,ansi}
HM450/<label><dir>.{bed,ansi}
manual/<lineage>/<region>[_<variant>].bed + shared region .ansi   # hand-defined
signature_master.tsv               # one row per whole-genome signature (auto + manual)
signature_master_hm450.tsv         # one row per HM450 label
LineageMarker.20260619.cm          # is_best-per-folder auto + all manual, fmt-s genome track
LineageMarker_20260619_hm450.rds   # SummarizedExperiment: HM450 probe beta x 65 cells
```

- `auto/<…>/*.bed` — the signature's selected CpG **locations** (the canonical,
  compact representation), in **beg0/end1** (0-based begin, 1-based end; a single
  CpG = `(C-1, C)`). A genome-wide binary `.cm` for any one signature is derived on
  demand:
  `awk 'NR==FNR{a[$1"_"$3]=1;next}{print (($1"_"($2+1)) in a)?1:0}' sig.bed <(zcat cpg_nocontig.bed.gz) | yame pack -f b - sig.cm`
- `auto/<…>/*.ansi` — colored decile inspection panel: the SEL track (labelled with
  the **full signature name**) then all 65 cells in lineage→sublineage order, target
  rows **underlined**, cropped to the selection ± 50 flanking CpGs.
- `HM450/<label>.bed` — the label's selected HM450 probes; `HM450/<label>.ansi` —
  real-beta (`L<0.34 M H>0.67`) panel, cells × probes, target rows underlined.
- `manual/` — see **Manual signatures** below.

## signature_master.tsv

Columns: `signature, cell_type, category (specific|broad|subgroup),
source (auto|manual), gene, chrm, region_beg, direction, n_cpg, is_best,
overlaps_manual, validated, matched_cellmarker`.

- `category`: `specific` = one cell type, `broad` = pan-lineage, `subgroup` =
  within-lineage subgroup.
- `is_best` flags the largest-`n_cpg` **auto** signature per folder; computed from
  `signature_master` alone, **independent of `LineageMarker` membership**.
- `LineageMarker.20260619.cm` holds one label per CpG: the `is_best` auto signature
  per folder (cell-type + pan + subgroup) **+ all manual** signatures, overlaps
  resolved **manual > auto** then master order, `MIN_CPG=5` re-enforced on the final
  track → **64 labels**. `overlaps_manual` names the manual signature an auto
  signature shares CpGs with (and loses to).

## Contents

- **auto/**: **57 folders** = cell types + `panEpith` / `panHema` / `panNeural` +
  **8 subgroup folders** (`subHema_{lymphoid,myeloid}`, `subNeural_{neuronal,glial}`,
  `subMesench_{skeletal,smoothfibro,endothelial}`, `subEpith_GI`). `panMesench`,
  `tongue_epi`, and `subEpith_nonGI` have no locus that survives the specificity
  refine and are dropped. **416 auto signatures** + **7 manual**.
- **`LineageMarker.20260619.cm`** → **64 labels** (is_best auto per folder + 7
  manual).
- **HM450**: **88 labels** (16 broad + 62 specific + 10 subgroup);
  `LineageMarker_20260619_hm450.rds` = **12,474 probe × 65 cell** β matrix. The
  broad `panEpith-` is the **hand-curated** set: the weak auto MRMP `panEpith-`
  (12 probes) is replaced by the 12 per-region manual signatures
  (e.g. `panEpith-.ELF3_chr1_202005002`, 41 probes total; `n_mrmp=0`), which take
  priority over auto. `panHema-`/`panNeural-` stay auto.
- `signature_master.tsv`: **423 rows** (416 auto + 7 manual).

## Manual signatures (`manual/`)

Hand-defined cell-type groups the automatic taxonomy cannot express. A manual
*region* may carry several bed **variants** (different selections at the same
locus): `<region>_<variant>.bed`. The variant is the single source of truth — it is
the SEL track label in the shared region `.ansi` and the `signature` name in
`signature_master.tsv`. Current set: `panEpith/` (ELF3, MIR200C, MIR200BA with
`_inclBasal`/`_exclBasal` variants), `panHema/` (PTPRC, WDFY4), `panNeural/` (OMG).
All are `source=manual`, `category=broad`, and **all** enter `LineageMarker`
(unlike auto, where only `is_best` does). Region matrices are extracted with
`yame rowsub -L <coords> sub.cg | yame unpack -a -f 5` (real β, NA if cov<5) — not
`hprint` — and rendered with one SEL track per bed. See the lab journal for the
regeneration commands.

## Methods (summary)

**Reference panel.** 65 cell types: 31 Loyfer per-cell-type pseudobulks
(`yame rowop -o musum`) + 34 Zhou "Major" pseudobulks, restricted to `include=1` in
the control sheet, across 5 lineages. 15 cross-dataset twins (e.g. NK_cell≡Hema_NK,
endothelial_cell≡Endo_Vsc) are collapsed by `equiv_group`. Adrenal cortex
(`Epi_AdrCtx`) is its own `Adrenal` lineage (it is methylated at pan-epithelial
markers, so it is not epithelial); endothelial cells are included (`Mesench`).

**MRMP discovery.** At each of 29,401,795 CpGs, keep sites with all 65 cell types at
depth ≥10, high-vs-low group mean Δ ≥0.6, not chrY; reduce to a 65-bit
methylated/unmethylated string; rank by genome-wide frequency, keep the top 10,000
(MRMPs). The MRMP set carries no lineage info — lineage/sublineage are applied only
downstream.

**Classification.** Each MRMP's minority ("auto") group, after twin-collapse, by
finest exact match: one cell type (**cell-type-specific**); its projection onto one
lineage equals one **sublineage** (**pan-sublineage** `sub<L>_<sub>`, other lineages
don't-care); ≥2 cell types all in one lineage (**pan-lineage** `pan{Epith,Hema,
Neural,Mesench}`); else **mixed**. Direction `-`/`+`.

**Curation + refine.** Per-label Win10k enrichment → top-5 merged windows → decile
base-curation; a lineage-free worst-case **specificity refine** removes
lineage-confounded CpGs. Subgroups curate against an explicit *rest-of-lineage*
outgroup (other lineages don't-care). Signatures with `< 5` CpGs are dropped.

**LineageMarker.** Per-signature CpG-location BEDs; the genome-wide `.cm` is built
once by concatenating the `is_best` auto + all manual beds, labelling by signature
name, resolving overlaps (manual > auto), and a single `yame pack -f s`. Panels
(`.ansi`) are rendered last (after gene/validation).

**HM450.** For each cell-type-specific / pan-lineage / pan-sublineage MRMP in the
top 1000 with ≥10 HM450 probes, continuous β (`yame unpack -f 1`) is taken. Broad
and subgroup labels are tightened by a per-probe **AUC rank-separation** filter
(target vs **all other cells**, AUC ≥0.99): a probe is kept only if its β cleanly
separates the target from every other cell — removing cross-lineage false positives
and heterogeneous targets. `panMesench` has no probe that separates all of
mesenchyme and is dropped. Cell-type labels are kept as-is.

## Pipeline

| step | script | output |
|---|---|---|
| 0 | `sa_00_reference.sh` | 65-cell `sub.cg` (Loyfer musum + Zhou Major) + cellorder |
| 1 | `sa_01_mrmp.sh` + `sa_01_def.R` | MRMP mask (Δmean ≥0.6, count==65, top-10000); def table |
| 2 | `sa_02_classify.R` | `mrmp_class.tsv` (cell-type / pan-sublineage / pan-lineage / mixed) |
| 3 | `sa_03_hm450_win.sh` | HM450 BED, Win10k set, HM450 MRMP probes |
| 4 | `sa_04_enrich_curate.sh` (+ `sa_refine_select.R`) | per-region base curation + specificity refine (cell-type/pan + subgroup passes) |
| 5 | `sa_05_auto.sh` | `auto/` beds + `signature_master.tsv` + `LineageMarker.<date>.cm` |
| 5h | `sa_05_hm450.sh` + `sa_05_hm450.R` | `HM450/` panels/beds + master + RDS |
| 6 | `sa_06_annotate_validate.R` | gene + CellMarker validation columns |
| 7 | `sa_render_panels.R` | `auto/` `.ansi` panels (rendered last) |
