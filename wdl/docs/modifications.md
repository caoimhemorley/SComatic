# SComatic modifications

The container used by these WDLs contains project-specific changes to [SComatic](https://github.com/cortes-ciriano-lab/SComatic):

- Fixed a Fisher's exact test bug.
- Extended `SingleCellGenotype.py` to report positional read-length information.
- Modified scripts to work with compressed input files.
- Added the custom `step2_pass_umi_library.py` for UMI-support assessment.
- Changed the single-cell genotype WDL default minimum base quality from 30 to 20, in line with BaseCount.

This is a functional summary, not a source-level diff. The source modifications are contained in the published `caoimhem/scomatic:0.3-amd64` image; the modified source tree is not supplied in this archive.
