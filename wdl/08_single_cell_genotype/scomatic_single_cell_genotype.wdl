version 1.0

workflow scomatic_singlecellgenotype {

  input {
    String sample_id

    File bam
    File bai
    File step2_pass

    File ref_fasta
    File ref_fai

    File metadata

    String alt_flag = "All"

    Int nprocs = 4
    Int bin = 50000
    Int min_bq = 20
    Int min_mq = 255
    Int max_dp = 8000
  }

  call single_cell_genotype {
    input:
      sample_id = sample_id,
      bam = bam,
      bai = bai,
      step2_pass = step2_pass,
      ref_fasta = ref_fasta,
      ref_fai = ref_fai,
      metadata = metadata,
      alt_flag = alt_flag,
      nprocs = nprocs,
      bin = bin,
      min_bq = min_bq,
      min_mq = min_mq,
      max_dp = max_dp
  }

  output {
    File genotype_calls = single_cell_genotype.genotype_calls
    File log_file = single_cell_genotype.log_file
  }
}

task single_cell_genotype {

  input {
    String sample_id

    File bam
    File bai
    File step2_pass

    File ref_fasta
    File ref_fai

    File metadata

    String alt_flag

    Int nprocs
    Int bin
    Int min_bq
    Int min_mq
    Int max_dp
  }

  command <<<
    set -euo pipefail

    mkdir -p input output logs tmp

    echo "Running SComatic SingleCellGenotype"
    echo "Sample: ~{sample_id}"

    ln -s "~{bam}" "input/~{sample_id}.bam"
    ln -s "~{bai}" "input/~{sample_id}.bam.bai"

    ln -s "~{ref_fasta}" input/genome.fa
    ln -s "~{ref_fai}" input/genome.fa.fai

    echo "Checking inputs..."

    test -s "input/~{sample_id}.bam"
    test -s "input/~{sample_id}.bam.bai"
    test -s input/genome.fa
    test -s input/genome.fa.fai
    test -s "~{metadata}"
    test -s "~{step2_pass}"

    echo "BAM:"
    ls -lh "input/~{sample_id}.bam"

    echo "Step2 PASS:"
    ls -lh "~{step2_pass}"

    echo "Metadata:"
    ls -lh "~{metadata}"

    if [[ "~{step2_pass}" == *.gz ]]; then
        gzip -dc "~{step2_pass}" > input/step2.pass.tsv
    else
        cp "~{step2_pass}" input/step2.pass.tsv
    fi

    test -s input/step2.pass.tsv

    echo "Starting SingleCellGenotype.py..."

    /opt/conda/envs/SComatic/bin/python \
      /opt/SComatic/scripts/SingleCellGenotype/SingleCellGenotype.py \
      --bam "input/~{sample_id}.bam" \
      --infile input/step2.pass.tsv \
      --ref input/genome.fa \
      --meta "~{metadata}" \
      --outfile "output/~{sample_id}.SingleCellGenotype.tsv" \
      --alt_flag "~{alt_flag}" \
      --nprocs ~{nprocs} \
      --bin ~{bin} \
      --min_bq ~{min_bq} \
      --min_mq ~{min_mq} \
      --max_dp ~{max_dp} \
      --tmp_dir tmp \
      2>&1 | tee "logs/~{sample_id}.SingleCellGenotype.log"

    test -s "output/~{sample_id}.SingleCellGenotype.tsv"

    echo "Checking output header..."
    head -1 "output/~{sample_id}.SingleCellGenotype.tsv"

    echo "Compressing output..."

    pigz -p ~{nprocs} \
      "output/~{sample_id}.SingleCellGenotype.tsv"

    test -s "output/~{sample_id}.SingleCellGenotype.tsv.gz"

    echo "SingleCellGenotype complete."

    ls -lh \
      "output/~{sample_id}.SingleCellGenotype.tsv.gz" \
      "logs/~{sample_id}.SingleCellGenotype.log"
  >>>

  output {
    File genotype_calls = "output/~{sample_id}.SingleCellGenotype.tsv.gz"
    File log_file = "logs/~{sample_id}.SingleCellGenotype.log"
  }

  runtime {
    docker: "caoimhem/scomatic:0.3-amd64"
    cpu: nprocs
    memory: "32 GB"
    disks: "local-disk 200 HDD"
    preemptible: 3
    bootDiskSizeGb: 30
  }
}
