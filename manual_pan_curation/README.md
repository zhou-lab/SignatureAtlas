# Manually-curated pan-lineage signatures (archived for comparison)

These 6 pan-lineage region signatures were **hand-curated** in earlier sessions and are **NOT**
produced by the reproducible pipeline (`labjournal/zhouw3/2026/tools/build_signature_atlas.sh`).
They are kept here only to compare against the auto-derived pan-lineage signatures that now live
in `whole_genome/panEpith/`, `panHema/`, `panNeural/`, `panMesench/`.

| folder | signature | gene |
|--------|-----------|------|
| panEpi   | ELF3, MIR200C, MIR200BA | epithelial |
| panLeuko | PTPRC, WDFY4            | leukocyte  |
| panNeu   | OMG                     | neural     |

The pipeline derives pan-lineage signatures objectively: MRMPs whose discordant cell types all
belong to one lineage (Epith / Hema / Neural / Mesench) are collapsed to that lineage's label and
their top Win10k-enriched regions curated, same as cell-type-specific markers. The manual set above
used a different, hand-picked gene-anchored approach; differences between the two are expected.
