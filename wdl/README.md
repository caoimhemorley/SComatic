# SComatic workflows

WDL workflows for identifying somatic single-nucleotide variants (SNVs) from single-nucleus RNA-seq data, based on [SComatic](https://github.com/cortes-ciriano-lab/SComatic).

This repository includes modifications to the original SComatic scripts, additional analysis functionality, and WDL workflows for running the pipeline. See [Modifications](docs/modifications.md) for a summary of changes.

## Requirements

- Aligned BAM files (`.bam`) and corresponding indexes (`.bai`).
- Cell barcode annotations mapping barcodes to cell types.
- Reference genome FASTA and index (GRCh38/hg38).
- A WDL execution engine, such as Cromwell, with access to the input files and Docker.
- The SComatic Docker image: [`caoimhem/scomatic:0.3-amd64`](https://hub.docker.com/r/caoimhem/scomatic).

The Docker image contains the SComatic scripts and required software dependencies, including the modifications used in these workflows.

## Main workflow

The main SComatic pipeline consists of five stages, run in the following order:

| Step | Workflow | Description |
|---|---|---|
| 1 | `01_splitbam` | Split BAM files by cell type. |
| 2 | `02_basecount` | Generate nucleotide counts for each cell type. |
| 3 | `03_merge_basecounts` | Merge base counts across cell types. |
| 4 | `04_basecellcalling_step1` | Identify candidate somatic variants. |
| 5 | `05_basecellcalling_step2` | Filter candidate variants and generate final variant calls. |

The workflows use the following cell-type categories: Astrocyte, Excitatory, Inhibitory, Microglia, Oligodendrocyte, OPC and Other.

Each workflow folder contains:

- The WDL file.
- A `wrapper.json` defining the workflow input mapping.
- An `inputs/example.csv` showing the expected input format.

The example CSVs use placeholder file paths and sample identifiers, which should be replaced with the appropriate dataset inputs.

For a description of the inputs and parameters for each workflow, see [WDL inputs and parameters](docs/parameters.md).

## Additional workflows

The following workflows provide additional analyses and can be run independently when their required inputs are available.

| Workflow | Description |
|---|---|
| `06_callable_sites` | Calculate callable sites. |
| `07_step2_pass_umi` | Assess UMI-level support for candidate variants. |
| `08_single_cell_genotype` | Genotype variants at single-cell resolution and obtain additional read-level information. |

## Docker

The workflows use the publicly available Docker image:

```bash
docker pull caoimhem/scomatic:0.3-amd64
```

Updated Docker images will be published to [Docker Hub](https://hub.docker.com/r/caoimhem/scomatic) with new version tags. Changes to scripts and workflows are tracked in this GitHub repository.

For details of the original SComatic method and software, refer to the [SComatic documentation](https://github.com/cortes-ciriano-lab/SComatic).