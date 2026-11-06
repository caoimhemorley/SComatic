version 1.0

workflow scomatic_basecount {
  input {
    String sample_id
    String cell_type
    File bam
    File bai
    File ref_fasta
    File ref_fai
    String chrom
    Int nprocs = 4
    Int bin_size = 50000
    Int min_dp = 5
    Int min_cc = 5
    Int min_bq = 20
    Int min_mq = 255
  }

  call base_cell_counter {
    input:
      sample_id = sample_id,
      cell_type = cell_type,
      bam = bam,
      bai = bai,
      ref_fasta = ref_fasta,
      ref_fai = ref_fai,
      chrom = chrom,
      nprocs = nprocs,
      bin_size = bin_size,
      min_dp = min_dp,
      min_cc = min_cc,
      min_bq = min_bq,
      min_mq = min_mq
  }

  output {
    File basecount_tsv_gz = base_cell_counter.basecount_tsv_gz
    File log_file = base_cell_counter.log_file
  }
}

task base_cell_counter {
  input {
    String sample_id
    String cell_type
    File bam
    File bai
    File ref_fasta
    File ref_fai
    String chrom
    Int nprocs
    Int bin_size
    Int min_dp
    Int min_cc
    Int min_bq
    Int min_mq
  }

  command <<<
    set -euo pipefail

    mkdir -p basecounts tmp logs

    /opt/conda/envs/SComatic/bin/python /opt/SComatic/scripts/BaseCellCounter/BaseCellCounter.py \
      --bam ~{bam} \
      --ref ~{ref_fasta} \
      --chrom ~{chrom} \
      --out_folder basecounts \
      --id ~{sample_id}.~{cell_type} \
      --nprocs ~{nprocs} \
      --bin ~{bin_size} \
      --min_dp ~{min_dp} \
      --min_cc ~{min_cc} \
      --min_bq ~{min_bq} \
      --min_mq ~{min_mq} \
      --tmp_dir tmp \
      > logs/~{sample_id}.~{cell_type}.BaseCellCounter.log 2>&1

    test -f "basecounts/~{sample_id}.~{cell_type}.tsv"
    pigz -p ~{nprocs} "basecounts/~{sample_id}.~{cell_type}.tsv"
  >>>

  output {
    File basecount_tsv_gz = "basecounts/~{sample_id}.~{cell_type}.tsv.gz"
    File log_file = "logs/~{sample_id}.~{cell_type}.BaseCellCounter.log"
  }

  runtime {
    docker: "caoimhem/scomatic:0.3-amd64"
    cpu: nprocs
    memory: "16 GB"
    disks: "local-disk 250 HDD"
    preemptible: 3
    bootDiskSizeGb: 30
  }
}
