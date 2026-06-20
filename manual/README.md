# Manually-curated signatures

These region signatures are **hand-defined** (cell-type groups chosen by hand) — **NOT** discovered by the automatic MRMP pipeline (`labjournal/zhouw3/2026/tools/sa_*`), though their selections are computed from the same 61-cell signal matrix. They capture groups the automatic lineage taxonomy cannot express, and (unlike auto signatures, where only the per-folder `is_best` enters) **all** of them are merged into `LineageMarker.20260619.cm` and indexed in `signature_master.tsv` (`source=manual`, `category=broad`).

## Region variants

A manual **region** (`<label>.<gene>_<chr>_<beg>`) may carry several bed **variants** — different selections at the same locus — named `<region>_<variant>.bed`. The variant suffix is the single source of truth: it is the SEL track label in the region's shared `.ansi` and the `signature` name in `signature_master.tsv`.

| folder | signatures | target group |
|--------|-----------|--------------|
| panEpith  | `ELF3`, `MIR200C` | pan-epithelial |
| panEpith  | `MIR200BA_inclBasal` / `MIR200BA_exclBasal` | epithelial **with** vs **without** breast-basal (myoepithelial) in the hypomethylated group |
| panHema   | `PTPRC`, `WDFY4` | leukocyte / hematopoietic |
| panNeural | `OMG` | neural |

**MIR200BA two-variant rationale.** At this locus breast basal (myoepithelial) cells do **not** follow the epithelial hypomethylation (smoothed decile ~5 vs ~0 for the other epithelial cell types). The two variants are computed directly from the signal and are **disjoint by construction**: a CpG enters `_inclBasal` if breast basal is *low* (part of the epithelial program) and `_exclBasal` if breast basal is *high* (excluded — "non-myoepithelial"). Hepatocyte and kidney are treated as **neutral** (no constraint) but still shown in the panel.

## Provenance / regeneration

Region matrices are extracted the proper way — `yame rowsub -R <cr> -L <coords> sub.cg` then `yame unpack -a -f 5` (real β, NA if cov<5), **not** `hprint` (visualization only) — via `tools/sa_extract_region.R`; the shared region panel (one SEL track per bed, label = bed filename) is rendered with `tools/sa_manual_panel.R`. The regeneration commands are in the lab journal: `labjournal/zhouw3/2026/20260619_signature_atlas.org` (section **Manual signatures**).
