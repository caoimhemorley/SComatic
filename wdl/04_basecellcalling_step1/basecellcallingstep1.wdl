version 1.0

workflow scomatic_step1 {
  input {
    String sample_id
    File merged_basecounts
    File ref_fasta

    Int min_cov = 5
    Int min_cells = 5
    Int min_ac_cells = 2
    Int min_ac_reads = 3
    Int max_cell_types = 4
    Int min_cell_types = 1

    Float fisher_cutoff = 0.001
    Float alpha1 = 0.260288007167716
    Float beta1 = 173.94711910763732
    Float alpha2 = 0.08354121346569514
    Float beta2 = 103.47683488327257

    Int nprocs = 4
  }

  call basecellcalling_step1 {
    input:
      sample_id = sample_id,
      merged_basecounts = merged_basecounts,
      ref_fasta = ref_fasta,
      min_cov = min_cov,
      min_cells = min_cells,
      min_ac_cells = min_ac_cells,
      min_ac_reads = min_ac_reads,
      max_cell_types = max_cell_types,
      min_cell_types = min_cell_types,
      fisher_cutoff = fisher_cutoff,
      alpha1 = alpha1,
      beta1 = beta1,
      alpha2 = alpha2,
      beta2 = beta2,
      nprocs = nprocs
  }

  output {
    File step1_calls = basecellcalling_step1.step1_calls
    File log_file = basecellcalling_step1.log_file
  }
}

task basecellcalling_step1 {
  input {
    String sample_id
    File merged_basecounts
    File ref_fasta

    Int min_cov
    Int min_cells
    Int min_ac_cells
    Int min_ac_reads
    Int max_cell_types
    Int min_cell_types

    Float fisher_cutoff
    Float alpha1
    Float beta1
    Float alpha2
    Float beta2

    Int nprocs
  }

  command <<<
    set -euo pipefail

    mkdir -p step1 logs

    test -s "~{merged_basecounts}"
    test -s "~{ref_fasta}"

    echo "Running SComatic BaseCellCalling Step 1 for ~{sample_id}"

    /opt/conda/envs/SComatic/bin/python \
      /opt/SComatic/scripts/BaseCellCalling/BaseCellCalling.step1.py \
      --infile "~{merged_basecounts}" \
      --outfile "step1/~{sample_id}" \
      --ref "~{ref_fasta}" \
      --min_cov ~{min_cov} \
      --min_cells ~{min_cells} \
      --min_ac_cells ~{min_ac_cells} \
      --min_ac_reads ~{min_ac_reads} \
      --max_cell_types ~{max_cell_types} \
      --min_cell_types ~{min_cell_types} \
      --fisher_cutoff ~{fisher_cutoff} \
      --alpha1 ~{alpha1} \
      --beta1 ~{beta1} \
      --alpha2 ~{alpha2} \
      --beta2 ~{beta2} \
      > "logs/~{sample_id}.BaseCellCalling.step1.log" 2>&1

    test -s "step1/~{sample_id}.calling.step1.tsv"

    pigz \
      -p ~{nprocs} \
      "step1/~{sample_id}.calling.step1.tsv"

    test -s "step1/~{sample_id}.calling.step1.tsv.gz"

    echo "Step 1 complete."
  >>>

  output {
    File step1_calls =
      "step1/${sample_id}.calling.step1.tsv.gz"

    File log_file =
      "logs/${sample_id}.BaseCellCalling.step1.log"
  }

  runtime {
    docker: "caoimhem/scomatic:0.3-amd64"
    cpu: nprocs
    memory: "32 GB"
    disks: "local-disk 500 HDD"
    preemptible: 3
    bootDiskSizeGb: 30
  }
}
