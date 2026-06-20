# SignatureAtlas

Whole-genome and HM450 DNA-methylation **marker signatures** for human cell types and lineages, derived from a 61-cell-type deep WGBS reference (30 Loyfer per-cell-type pseudobulks + 31 Zhou "Major" pseudobulks). Markers are Most-Recurrent Methylation Patterns (MRMPs) curated into compact genomic regions, labelled `<cell_type>` (cell-type-specific) or `pan<Lineage>` (pan-lineage), with direction `-` = hypomethylated and `+` = hypermethylated in the target cell type(s).

## Layout

```
whole_genome/<cell_type|panLineage>/<label><dir>.<gene>_<chr>_<beg>.{bed,ansi}
HM450/<label><dir>.{bed,ansi}
signature_master.tsv              # one row per whole-genome signature (auto + manual)
signature_master_hm450.tsv        # one row per HM450 label
LineageMarker.20260619.cm         # is_best-per-folder auto + all manual signatures, fmt-s genome track
LineageMarker_20260619_hm450.rds  # SummarizedExperiment: HM450 probe beta x 61 cells
manual/                           # hand-defined signatures (NOT pipeline output; source=manual)
```

- `whole_genome/<…>/*.bed` — the signature's selected CpG **locations** (the canonical, compact representation), in **beg0/end1** (0-based begin, 1-based end; a single CpG = `(C-1, C)`). A genome-wide binary `.cm` for any one signature is derived on demand:
  `awk 'NR==FNR{a[$1"_"$3]=1;next}{print (($1"_"($2+1)) in a)?1:0}' sig.bed <(zcat cpg_nocontig.bed.gz) | yame pack -f b - sig.cm`
- `whole_genome/<…>/*.ansi` — colored decile inspection panel: the SEL track (labeled with the **full signature name**) then all 61 cells in lineage order, cropped to the selection ± 50 flanking CpGs.
- `HM450/<label>.bed` — all HM450 probes mapping to that label's MRMPs; `HM450/<label>.ansi` — real-beta (`L<0.34 M H>0.67`) panel, cells × probes.
- `manual/` — see **Manual signatures** below.

`signature_master.tsv` columns: `signature, cell_type, category (specific|broad), source (auto|manual), gene, chrm, region_beg, direction, n_cpg, is_best, overlaps_manual, validated, matched_cellmarker`.
- `is_best` flags the largest-`n_cpg` **auto** signature per folder; it is computed from `signature_master` alone and is **independent of `LineageMarker` membership**.
- `LineageMarker` holds one label per CpG with **manual priority** (manual wins every shared CpG), then `MIN_CPG=5` is re-enforced on the final track. `overlaps_manual` names the manual signature an auto signature shares CpGs with. An `is_best` auto signature whose CpGs are mostly claimed by a manual signature drops below 5 CpGs here and is **removed from `LineageMarker`** — this happens to `panEpith-.AL390719.2_chr1_1100000` (a duplicate-discovery of the manual MIR200BA locus, left with 1 CpG). So `LineageMarker` has **52 labels** = 46 `is_best` auto − 1 + 7 manual.

## Contents

- **46 folders** = 43 cell types + panEpith / panHema / panNeural. (`panMesench` and `tongue_epi` have no locus that survives the specificity refine and are dropped from the auto set.)
- **355 auto whole-genome signatures** + **7 manual** signatures (`manual/`). Signatures with `< 5` selected CpGs after refine are dropped.
- **`LineageMarker.20260619.cm`** = the `is_best` auto signature per folder + **all 7 manual** signatures (manual always enters; auto only via `is_best`), conflict-resolved (manual wins) and `MIN_CPG=5`-filtered on the final track to **52 labels**, each keyed by its full signature name.
- **67 HM450 labels**; RDS = 17,457 probe × 61 cell beta matrix.
- representative genes validated against CellMarker 2.0 (`validated` column).

## Manual signatures (`manual/`)

Hand-defined cell-type groups the automatic lineage taxonomy cannot express. A manual *region* may carry several bed **variants** (different selections at the same locus): `<region>_<variant>.bed`. The variant is the single source of truth — it is the SEL track label in the shared region `.ansi` and the `signature` name in `signature_master.tsv`. Current set:

| folder | signatures | target group |
|--------|-----------|--------------|
| panEpith | `ELF3`, `MIR200C` | pan-epithelial |
| panEpith | `MIR200BA_inclBasal` / `MIR200BA_exclBasal` | epithelial *with* / *without* breast-basal (myoepithelial) in the hypomethylated group — two **disjoint** variants (breast basal low vs high); hepatocyte & kidney are neutral (shown, not constrained) |
| panHema | `PTPRC`, `WDFY4` | leukocyte / hematopoietic |
| panNeural | `OMG` | neural |

All are `source=manual`, `category=broad`, and **all** enter `LineageMarker` (unlike auto). Region matrices for manual panels are extracted with `yame rowsub -L <coords> sub.cg | yame unpack -a -f 5` (real β, NA if cov<5) — not `hprint` — and rendered with one SEL track per bed (label = bed filename). See the lab journal for the regeneration commands.

## Methods

**Reference panel.** 61 cell types from two public sources: 30 Loyfer per-cell-type pseudobulks (`yame rowop -o musum`) + 31 Zhou "Major" pseudobulks, restricted to `include=1` in the control sheet. 14 cross-dataset twins (e.g. NK_cell≡Hema_NK) are collapsed by `equiv_group` into one folder each.

**MRMP discovery.** At each of 29,401,795 CpGs, keep sites with all 61 cell types at depth ≥10, high-vs-low group mean Δ ≥0.6, not chrY; reduce to a 61-bit methylated/unmethylated string (`yame rowop -o stat`/`binstring`); rank by genome-wide frequency, keep the top 10,000 (MRMPs).

**Classification.** Each MRMP's minority ("auto") group, after twin-collapse, resolves to one cell type (**cell-type-specific**), ≥2 cell types in one lineage (**pan-lineage** `pan{Epith,Hema,Neural,Mesench}`), or **mixed** (not curated). Direction `-`/`+`.

**Region enrichment + base curation.** MRMPs are relabelled by class and enriched against fixed 10 kb windows (`yame summary`, Log2 odds, KYCG); per label, windows overlapping ≥5 of its MRMP CpGs are kept, merged, top 5 curated. `curate_window.sh` extracts each region's decile matrix once (`yame hprint -g`) and base-selects a CpG if it is unmethylated (β<0.4) in ≥90% of the target hypo group **and** methylated (β≥0.6) in ≥80% of the complementary hyper group (±100-CpG flank).

**Specificity refine.** The pooled base test is diluted by panel composition — a CpG hypomethylated across the target's *entire lineage* can still clear the 80% bar (e.g. a B-cell marker that is really pan-leukocyte). A lineage-free filter (`sa_refine_select.R`) smooths each cell's per-CpG decile signal over a centered 5-CpG window and keeps a selected CpG only if the **worst-case** outgroup cell is still separated (hypo: min smoothed ≥5; hyper: max smoothed ≤4). Subtractive; signatures with `<5` CpGs are dropped. (`panMesench`, `tongue_epi` do not survive.)

**Whole-genome output.** Each signature is stored as its CpG-location BED (beg0/end1); the genome-wide `.cm` is derived on demand. `LineageMarker.<date>.cm` (one label per CpG) is built by concatenating the per-folder `is_best` auto beds + all manual beds, labeling each by signature name, resolving overlaps (manual > auto), and a single `yame pack -f s` — turning ~356 genome scans into one. Panels (`.ansi`) are rendered **last** (after gene/validation) so each SEL track carries the full signature name.

**HM450.** For each cell-type-specific / pan-lineage MRMP in the top 1000 with ≥10 HM450 probes, all probes are kept; continuous β (`yame unpack -f 1`) is written per label and shipped as a SummarizedExperiment (17,457 probe × 61 cell).

**Validation.** Each signature's representative gene is cross-checked against CellMarker 2.0; supporting evidence (`validated`/`matched_cellmarker`), not a filter.

## Reproducibility

Regenerated from raw `.cg` by self-contained staged scripts in the lab journal (`labjournal/zhouw3/2026/tools/`). There is no end-to-end driver — run the steps in order (each is one hop; see the org for the exact commands and the env preamble):

| step | script | output |
|------|--------|--------|
| setup | `sa_setup_cellmarker.R` (curl in the org) | CellMarker 2.0 → references annotation dir |
| 0 | `sa_00_reference.sh` | 61-cell-type `sub.cg` (Loyfer musum + Zhou Major) |
| 1 | `sa_01_mrmp.sh` + `sa_01_def.R` | MRMP mask (Δmean ≥0.6, count==61, top-10000); def table |
| 2 | `sa_02_classify.R` | cell_type_specific / pan_lineage / mixed per MRMP |
| 3 | `sa_03_hm450_win.sh` | HM450 BED, Win10k set, HM450 MRMP probes |
| 4 | `sa_04_enrich_curate.sh` | per-label enrichment → base-curated region selections + decile matrices |
| 4b | `sa_refine_select.R` | specificity refine (filter only) of each region's selection |
| 5 | `sa_05_wholegenome.sh` | `whole_genome/*.bed`, `signature_master.tsv`, `LineageMarker.cm` |
| 5h | `sa_05_hm450.sh` + `sa_05_hm450.R` | HM450 real-β probe matrix; `HM450/`, `signature_master_hm450.tsv`, RDS |
| 6 | `sa_06_annotate_validate.R` | gene + CellMarker validation columns |
| 7 | `sa_render_panels.R` | `whole_genome/*.ansi` (rendered last; SEL track = signature name) |

Inputs: raw per-sample `.cg` under `/mnt/isilon/zhou_lab/projects/20230727_all_public_WGBS/hg38`, the control master sheet `labmeta/tsv/20260619_signature_atlas_celltypes.tsv`, gene coordinates (gencode.v36), the HM450 manifest, and CellMarker 2.0. See `labjournal/zhouw3/2026/20260619_signature_atlas.org` for the full record.
