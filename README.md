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
LineageMarker.20260619.cm     # best auto signature/folder + manual signatures, fmt-s genome track
LineageMarker_20260619_hm450.rds  # SummarizedExperiment: HM450 probe beta x 61 cells
marker_gene_validation.tsv    # CellMarker cross-check of representative genes
manual/                       # hand-curated signatures (NOT pipeline output; source=manual)
```

- `whole_genome/<…>/*.cm` — packed fmt-b genome-wide selection track (1 = selected CpG).
- `whole_genome/<…>/*.ansi` — colored decile panel (SEL row, target cell type(s), then all 61
  cells in lineage order) for visual inspection.
- `HM450/<label>.bed` — all HM450 probes mapping to that label's MRMPs.
- `HM450/<label>.ansi` — real-beta (`L<0.34 M H>0.67`) panel, cells × probes.
- `manual/` — hand-curated signatures kept aside; they use the same naming as auto signatures, are
  merged into `LineageMarker.20260619.cm`, and are indexed in `signature_master.tsv` (`source=manual`).

`signature_master.tsv` columns: `signature, cell_type, category (specific|broad),
source (auto|manual), gene, chrm, region_beg, direction, n_cpg, is_best, overlaps_manual,
validated, matched_cellmarker`.
- `is_best` flags the largest-`n_cpg` **auto** signature per folder. It is computed from
  `signature_master` alone and is **independent of `LineageMarker` membership**.
- `overlaps_manual` names the manual signature an auto signature shares CpGs with. `LineageMarker`
  holds one label per CpG with **manual priority**, so manual wins every shared CpG. If an auto
  signature's CpGs are **entirely** claimed by a manual signature, it contributes nothing and is
  **absent from `LineageMarker` even when it is `is_best`** — e.g. `panEpith-.AL390719.2_chr1_1100000`
  (all 59 CpGs ceded to `panEpith-.MIR200BA_chr1_1160292`; still `is_best=1` in the master). This is
  why `LineageMarker` has 53 keys rather than 48 best-auto + 6 manual = 54.

## Contents

- **48 folders** = 44 cell types + 4 pan-lineages (panEpith, panHema, panNeural, panMesench).
- **381 auto whole-genome signatures** (362 cell-type-specific, 19 pan-lineage) + **6 manual**
  signatures (`manual/`). `is_best` flags the highest-CpG auto signature per folder. (Signatures
  with 0 selected CpGs are filtered out at generation.)
- **`LineageMarker.20260619.cm`** = the best auto signature per folder + the 6 manual signatures,
  conflict-resolved (manual wins) to **53 keys**, each keyed by its full signature name.
- **67 HM450 labels**; RDS = 17,457 probe × 61 cell beta matrix.
- representative genes validated against CellMarker 2.0 (`validated` column).

## Methods

**Reference methylome panel.** A deep WGBS reference of 61 human cell types is assembled from two
public sources. For the 30 Loyfer *et al.* cell types, all samples of each cell type are pooled into
a single high-depth pseudobulk (`yame rowop -o musum`); for the 31 Zhou "Major" cell types the
existing pseudobulks are used as-is. The two are concatenated and restricted to the 61 cell types
flagged `include=1` in the control master sheet. Fourteen cell types measured in both datasets
(e.g. NK_cell≡Hema_NK, oligodendrocyte≡Glia_Oligo, pancreas_acinar≡Epi_Aci) are declared as
`equiv_group` twins and treated as a single biological entity throughout — removing spurious calls
between technical replicates of the same cell type and yielding better-ranked markers.

**Binary methylation patterns and MRMP discovery.** At each of the 29,401,795 genomic CpGs the
per-cell-type methylation level is summarised (`yame rowop -o stat`) and the site retained only if
(i) all 61 cell types have depth ≥10 (M+U), (ii) the high- vs low-group mean difference is ≥0.6,
and (iii) it is not on chrY. Each retained CpG is reduced to a 61-bit string of methylated vs
unmethylated cell types (`yame rowop -o binstring`). These **Most-Recurrent Methylation Patterns
(MRMPs)** are ranked by genome-wide frequency, and the 10,000 most frequent are kept — so signatures
reflect methylation programs shared by many loci rather than isolated single-CpG events.

**Signature classification.** Each MRMP's minority ("auto") group defines its discordant cell types.
After collapsing `equiv_group` twins, an MRMP resolving to one cell type is **cell-type-specific**;
one with ≥2 cell types all in one lineage (Epith / Hema / Neural / Mesench) is **pan-lineage**
(`pan{Epith,Hema,Neural,Mesench}`); the rest are **mixed** and not curated. Direction is `-` (target
hypomethylated) or `+` (hypermethylated). Pan-lineage signatures are thus derived objectively, in the
same framework as cell-type markers.

**Region enrichment and curation.** MRMPs are relabelled by class label and tested for enrichment
against fixed 10 kb windows (`yame summary`, Log2 odds ratio; KYCG framework). Per label, windows
overlapping ≥5 of its MRMP CpGs are kept, adjacent windows merged, and the top 5 regions curated.
Within a region a CpG is selected (`curate_window.sh`, via `yame hprint -g` deciles) if it is
unmethylated (β<0.4) in ≥90% of the target's hypo group **and** methylated (β≥0.6) in ≥80% of the
complementary hyper group; a ±100-CpG flank avoids clipping boundary DMRs.

**Whole-genome assembly.** Each curated CpG set is packed into a genome-wide track (`.cm`) and
rendered as a colored decile panel (`.ansi`), annotated with its locally dominant gene (gencode.v36,
±10 kb). The largest signature per folder is flagged `is_best`; their union forms
`LineageMarker.<date>.cm`.

**HM450 signatures.** For each cell-type-specific / pan-lineage MRMP in the top 1000 with ≥10 HM450
probes, all probes are retained; continuous β across the 61 cell types (`yame unpack -f 1`) is
written per label as a real-β panel + probe BED, and shipped whole as a SummarizedExperiment
(`LineageMarker_<date>_hm450.rds`, 17,457 probe × 61 cell).

**Marker-gene validation.** Each signature's representative gene is cross-checked against CellMarker
2.0 (cell-type→marker-name regex); results populate the `validated` / `matched_cellmarker` columns.
Because a methylation DMR need not coincide with the canonical expression-marker gene, validation is
supporting evidence, not a filter.

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
