# SignatureAtlas

Whole-genome and HM450 DNA-methylation **marker signatures** for human cell types and lineages,
derived from a 61-cell-type deep WGBS reference (30 Loyfer per-cell-type pseudobulks + 31 Zhou
"Major" pseudobulks). Markers are Most-Recurrent Methylation Patterns (MRMPs) curated into compact
genomic regions, labelled `<cell_type>` (cell-type-specific) or `pan<Lineage>` (pan-lineage), with
direction `-` = hypomethylated and `+` = hypermethylated in the target cell type(s).

## Layout

```
whole_genome/<cell_type|panLineage>/<label><dir>.<gene>_<chr>_<beg>.{cm,ansi}
HM450/<label><dir>.{bed,ansi}
signature_master.tsv          # one row per whole-genome signature
signature_master_hm450.tsv    # one row per HM450 label
LineageMarker.20260619.cm     # best signature per label, packed fmt-s genome-wide track
LineageMarker_20260619_hm450.rds  # SummarizedExperiment: HM450 probe beta x 61 cells
marker_gene_validation.tsv    # CellMarker cross-check of representative genes
manual_pan_curation/          # archived hand-curated pan sets (NOT pipeline output)
```

- `whole_genome/<…>/*.cm` — packed fmt-b genome-wide selection track (1 = selected CpG).
- `whole_genome/<…>/*.ansi` — colored decile panel (SEL row, target cell type(s), then all 61
  cells in lineage order) for visual inspection.
- `HM450/<label>.bed` — all HM450 probes mapping to that label's MRMPs.
- `HM450/<label>.ansi` — real-beta (`L<0.34 M H>0.67`) panel, cells × probes.

`signature_master.tsv` columns: `signature, cell_type, category (specific|broad), gene, chrm,
region_beg, direction, n_cpg, n_hm450, is_best, validated, matched_cellmarker`.

## Contents

- **48 folders** = 44 cell types + 4 pan-lineages (panEpith, panHema, panNeural, panMesench).
- **415 whole-genome signatures** (375 cell-type-specific, 40 pan-lineage); `is_best` flags the
  highest-CpG signature per folder (49 best, packed into `LineageMarker.20260619.cm`).
- **67 HM450 labels**; RDS = 17,457 probe × 61 cell beta matrix.
- 78 signatures' representative genes are validated against CellMarker 2.0.

## Reproducibility

The entire atlas is regenerated from raw `.cg` by a staged pipeline in the lab journal
(`labjournal/zhouw3/2026/tools/`), driven by `build_signature_atlas.sh`:

| stage | script | output |
|-------|--------|--------|
| setup | `sa_setup_cellmarker.sh` | CellMarker 2.0 → references annotation dir |
| 0 | `sa_00_reference.sh` | 61-cell-type `sub.cg` (Loyfer musum + Zhou Major) |
| 1 | `sa_01_mrmp.sh` | MRMP mask + def (delta_mean ≥0.6, count==61, top-10000) |
| 2 | `sa_02_classify.sh` | cell_type_specific / pan_lineage / mixed per MRMP |
| 3 | `sa_03_hm450_win.sh` | HM450 BED, Win10k set, HM450 MRMP probes |
| 4 | `sa_04_enrich_curate.sh` | per-label Win10k enrichment → curated region signatures |
| 5 | `sa_05_wholegenome.sh` | `whole_genome/`, masters, `LineageMarker.cm` |
| 5 | `sa_05_hm450.sh` | `HM450/`, `signature_master_hm450.tsv`, RDS (real beta) |
| 6 | `sa_06_annotate_validate.sh` | gene + CellMarker validation columns |

Inputs: raw per-sample `.cg` under `/mnt/isilon/zhou_lab/projects/20230727_all_public_WGBS/hg38`,
the control master sheet `labmeta/tsv/20260619_signature_atlas_celltypes.tsv` (cell-type →
dataset/lineage/equiv_group/include mapping), gene coordinates (gencode.v36), HM450 manifest,
and CellMarker 2.0. See `labjournal/zhouw3/2026/20260619_signature_atlas.org` for the full record.
