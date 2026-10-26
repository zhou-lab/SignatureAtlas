# SignatureAtlas releases

Each release is a git tag; the data files carry no version. Counts are for
the files at that tag: contrasts in `defs.tsv`, signed signatures, whole-genome
regions and CpGs in `wg.bed`, and HM450 probes and signatures in `hm450.bed`.

## v7 (2026-10-25)

119 contrasts · 205 signatures · 1,254 regions / 48,842 CpGs · 11,885 HM450
probes / 141 signatures

- Two blood references join the tissue atlas. Blood 5modC adds 10 contrasts
  (`*_Blood`; 95 regions, 5,865 CpGs) from sorted immune cells (EM-seq) and
  BLUEPRINT. Blood 5hmC adds 7 contrasts (`*_hmC`; 67 regions, 11,650 CpGs)
  from ACE-seq.
- `defs.tsv` has one row per contrast. In-group and out-group are lists of
  sample IDs (`<dataset>:<sample>`). New columns give the mark (5modC or
  5hmC) and the calling thresholds.
- `samples.tsv` (sample to cell type) and `datasets.tsv` (one row per source
  dataset) replace `reference_panel_samples.tsv`.
- The tissue regions, CpGs, HM450 probes and HM450 panels are unchanged from
  v6. CpGs are unique within each reference; 622 appear in two references.
- Rebuilds now break an exact margin tie between signatures by the MRMP's
  rank in that build, so renaming a contrast never moves a CpG. The blood
  references used this rule.

## v6 (2026-10-07)

102 contrasts · 171 signatures · 1,092 regions / 31,327 CpGs · 11,885 HM450
probes / 141 signatures

- Four myeloid contrasts renamed; nothing else changed:
  `Myeloid_MGC_Mast` → `Pan_Myeloid`, `Myeloid` → `Mature_Myeloid`,
  `Myeloid_Mast` → `Mature_Myeloid_Mast`, `Myeloid_Microglia` →
  `Mature_Myeloid_Microglia`.
- Adds `reference_panel_samples.tsv`, the raw samples behind each of the 87
  pseudobulks.
- First release with the browser at <https://zhou-lab.github.io/SignatureAtlas/>.

## v5 (2026-06-29)

102 contrasts · 171 signatures · 1,092 regions / 31,327 CpGs · 11,885 HM450
probes / 141 signatures

- Refines the myeloid contrasts. Mast cells are methylation-divergent and
  shrank the myeloid HM450 projection. A mast-free `Myeloid` (7 → 47 probes)
  and `Myeloid_Microglia` (19 → 78) are added. The mast-containing sets are
  kept as `Myeloid_Mast` and `Myeloid_MGC_Mast`.

## v4 (2026-06-29)

100 contrasts · 169 signatures · 1,075 regions / 31,122 CpGs · 11,760 HM450
probes / 139 signatures

- Adds 12 organ and lineage contrasts found among the top unannotated MRMPs:
  `Neural`, `Neural_CNS`, `Pancreas`, `Pancreatic_Exocrine`,
  `Respiratory_Thyroid`, `Hepato_Acinar`, `Mullerian_Renal`, `Mullerian_Epi`,
  `Striated_Muscle`, `Nerve_Fibroblast`, `Urogenital_Epi`, `Myelinating_Glia`.

## v3 (2026-06-29)

88 contrasts · 149 signatures · 957 regions / 28,306 CpGs · 11,166 HM450
probes / 120 signatures

- Regions come from direct clustering of each signature's MRMP CpGs, plus
  the HM450-probe-dense loci seeded on their probe span.
- Best-separation dedup: a CpG claimed by several signatures stays in the
  one with the widest margin, so no coordinate appears twice.
- Files lose the `v2_` prefix; the release is the git tag.

## v2 (2026-06-28)

83 contrasts · 144 signatures · 582 regions / 24,782 CpGs · 11,105 HM450
probes / 112 signatures

- First 87-cell-type reference: Loyfer, Zhou, BLUEPRINT and Tian
  pseudobulks. Signatures named by `signature_id` contrasts. Regions chosen
  by 10 kb window enrichment. Files carry a `v2_` prefix.
