# Manually-curated pan-lineage signatures (archived for comparison)

These 6 pan-lineage region signatures were **hand-curated** in earlier sessions and are **NOT**
produced by the reproducible pipeline (`labjournal/zhouw3/2026/tools/build_signature_atlas.sh`).
They are kept here to compare against the auto-derived pan-lineage signatures in
`whole_genome/panEpith/`, `panHema/`, `panNeural/`, `panMesench/`.

Folder/signature names use the **same lineage nomenclature as the auto signatures** (`panEpith`,
`panHema`, `panNeural`) so the two sets are directly comparable:

| folder | signatures | lineage |
|--------|-----------|---------|
| panEpith  | ELF3, MIR200C, MIR200BA | epithelial |
| panHema   | PTPRC, WDFY4            | leukocyte / hematopoietic |
| panNeural | OMG                     | neural |

These 6 are **also merged into `LineageMarker.20260619.cm`** and indexed in `signature_master.tsv`
(`category = manual_pan`). They take **priority over auto signatures** where CpGs overlap; each
affected auto signature records the manual signature it overlaps in the `overlaps_manual` column.

The pipeline derives pan-lineage signatures objectively: MRMPs whose discordant cell types all
belong to one lineage (Epith / Hema / Neural / Mesench) are collapsed to that lineage's label and
their top Win10k-enriched regions curated, same as cell-type-specific markers. The manual set above
used a different, hand-picked gene-anchored approach; differences between the two are expected.
