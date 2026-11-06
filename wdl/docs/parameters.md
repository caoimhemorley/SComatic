# WDL inputs and parameters

The WDLs reference the public image `caoimhem/scomatic:0.3-amd64`. Values below are the defaults in the packaged WDLs. **Required** means the workflow input has no WDL default. An input supplied at submission can override a default. The input CSVs show batch columns; options absent from those CSVs use WDL defaults unless provided separately.

The main workflow consists of stages 1–5. Callable sites, UMI support and single-cell genotyping are separate workflows, with no required order between them.

## 01_splitbam — `SplitBamWorkflow`

| Input | Type | WDL default |
|---|---|---|
| `bam` | `File` | Required |
| `bai` | `File` | Required |
| `metadata_tsv` | `File` | Required |
| `sample_id` | `String` | Required |
| `min_MQ` | `Int` | `255` |
| `max_nM` | `Int` | `5` |
| `max_NH` | `Int` | `1` |
| `n_trim` | `Int` | `0` |
| `cpu` | `Int` | `4` |
| `memory_gb` | `Int` | `16` |

**Outputs:** Cell-type BAMs and indexes (seven cell types), split report.

## 02_basecount — `scomatic_basecount`

| Input | Type | WDL default |
|---|---|---|
| `sample_id` | `String` | Required |
| `cell_type` | `String` | Required |
| `bam` | `File` | Required |
| `bai` | `File` | Required |
| `ref_fasta` | `File` | Required |
| `ref_fai` | `File` | Required |
| `chrom` | `String` | Required |
| `nprocs` | `Int` | `4` |
| `bin_size` | `Int` | `50000` |
| `min_dp` | `Int` | `5` |
| `min_cc` | `Int` | `5` |
| `min_bq` | `Int` | `20` |
| `min_mq` | `Int` | `255` |

**Outputs:** Per-cell-type compressed base-count TSV, log.

## 03_merge_basecounts — `scomatic_mergebasecounts_v2`

| Input | Type | WDL default |
|---|---|---|
| `sample_id` | `String` | Required |
| `astrocyte` | `File` | Required |
| `excitatory` | `File` | Required |
| `inhibitory` | `File` | Required |
| `microglia` | `File` | Required |
| `oligodendrocyte` | `File` | Required |
| `opc` | `File` | Required |
| `other` | `File` | Required |
| `nprocs` | `Int` | `4` |

**Outputs:** Merged compressed base-count TSV, log.

## 04_basecellcalling_step1 — `scomatic_step1`

| Input | Type | WDL default |
|---|---|---|
| `sample_id` | `String` | Required |
| `merged_basecounts` | `File` | Required |
| `ref_fasta` | `File` | Required |
| `min_cov` | `Int` | `5` |
| `min_cells` | `Int` | `5` |
| `min_ac_cells` | `Int` | `2` |
| `min_ac_reads` | `Int` | `3` |
| `max_cell_types` | `Int` | `4` |
| `min_cell_types` | `Int` | `1` |
| `fisher_cutoff` | `Float` | `0.001` |
| `alpha1` | `Float` | `0.260288007167716` |
| `beta1` | `Float` | `173.94711910763732` |
| `alpha2` | `Float` | `0.08354121346569514` |
| `beta2` | `Float` | `103.47683488327257` |
| `nprocs` | `Int` | `4` |

**Outputs:** Step 1 calls, log.

## 05_basecellcalling_step2 — `scomatic_step2`

| Input | Type | WDL default |
|---|---|---|
| `sample_id` | `String` | Required |
| `step1_calls` | `File` | Required |
| `min_distance` | `Int` | `5` |
| `nprocs` | `Int` | `2` |

**Outputs:** Step 2 calls and PASS-only calls, log.

## 06_callable_sites — `scomatic_callable_sites`

| Input | Type | WDL default |
|---|---|---|
| `sample_id` | `String` | Required |
| `step1_calls` | `File` | Required |
| `max_cov` | `Int` | `150` |
| `min_cell_types` | `Int` | `1` |

**Outputs:** Callable sites output, log.

## 07_step2_pass_umi — `Step2PassUMI`

| Input | Type | WDL default |
|---|---|---|
| `sample_id` | `String` | Required |
| `step2_calls` | `File` | Required |
| `astrocyte_bam` | `File` | Required |
| `astrocyte_bai` | `File` | Required |
| `excitatory_bam` | `File` | Required |
| `excitatory_bai` | `File` | Required |
| `inhibitory_bam` | `File` | Required |
| `inhibitory_bai` | `File` | Required |
| `microglia_bam` | `File` | Required |
| `microglia_bai` | `File` | Required |
| `oligodendrocyte_bam` | `File` | Required |
| `oligodendrocyte_bai` | `File` | Required |
| `opc_bam` | `File` | Required |
| `opc_bai` | `File` | Required |
| `other_bam` | `File` | Required |
| `other_bai` | `File` | Required |
| `max_records` | `Int` | `1000000000` |
| `min_bq` | `Int` | `20` |
| `min_mq` | `Int` | `255` |
| `max_depth` | `Int` | `8000` |

**Outputs:** Cell-type and variant UMI summaries, molecule evidence, log.

## 08_single_cell_genotype — `scomatic_singlecellgenotype`

| Input | Type | WDL default |
|---|---|---|
| `sample_id` | `String` | Required |
| `bam` | `File` | Required |
| `bai` | `File` | Required |
| `step2_pass` | `File` | Required |
| `ref_fasta` | `File` | Required |
| `ref_fai` | `File` | Required |
| `metadata` | `File` | Required |
| `alt_flag` | `String` | `"All"` |
| `nprocs` | `Int` | `4` |
| `bin` | `Int` | `50000` |
| `min_bq` | `Int` | `20` |
| `min_mq` | `Int` | `255` |
| `max_dp` | `Int` | `8000` |

**Outputs:** Compressed single-cell genotype TSV, log.

## Notes on settings

- `min_bq` is minimum Phred base quality; `min_mq`/`min_MQ` is minimum mapping quality. The BaseCount and UMI WDLs use Q20, and the SingleCellGenotype WDL in this package also defaults to Q20.
- `min_dp`/`min_cov` and `min_cc`/`min_cells` concern read depth and cell coverage/support, depending on the stage; check the relevant script for the exact computation.
- The alpha/beta values in Step 1 are beta-binomial model parameters. Their defaults are retained from the supplied WDL.
- Step 2 uses the hg38 RNA-editing database and PoN locations specified within the WDL; these are not exposed as workflow input parameters.
- Step2PassUMI relies on the project-specific script. Its molecule counts should not be assumed to equal read counts or independent genotyping calls.
- CPU, memory, disk and preemption runtime settings are declared inside the WDL task `runtime` blocks; the table above covers workflow-level inputs, not every execution-backend setting.
