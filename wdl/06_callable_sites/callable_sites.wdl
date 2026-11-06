version 1.0

workflow scomatic_callable_sites {

  input {
    String sample_id
    File step1_calls

    Int max_cov = 150
    Int min_cell_types = 1
  }

  call get_callable_sites {
    input:
      sample_id = sample_id,
      step1_calls = step1_calls,
      max_cov = max_cov,
      min_cell_types = min_cell_types
  }

  output {
    File callable_sites = get_callable_sites.callable_sites
    File log_file = get_callable_sites.log_file
  }
}

task get_callable_sites {

  input {
    String sample_id
    File step1_calls

    Int max_cov
    Int min_cell_types
  }

  command <<<
    set -euo pipefail

    mkdir -p output logs

    test -s "~{step1_calls}"

    echo "Running GetAllCallableSites for ~{sample_id}"
    echo "Compressed input: ~{step1_calls}"
    echo "max_cov: ~{max_cov}"
    echo "min_cell_types: ~{min_cell_types}"

    echo "Decompressing Step1 input..."
    gzip -dc "~{step1_calls}" > step1_input.tsv

    test -s step1_input.tsv

    echo "Decompressed input:"
    ls -lh step1_input.tsv

    /opt/conda/envs/SComatic/bin/python \
      /opt/SComatic/scripts/GetCallableSites/GetAllCallableSites.py \
      --infile step1_input.tsv \
      --outfile "output/~{sample_id}" \
      --max_cov ~{max_cov} \
      --min_cell_types ~{min_cell_types} \
      > "logs/~{sample_id}.callable_sites.log" 2>&1

    echo "Generated files:"
    ls -lh output/

    test "$(find output -maxdepth 1 -type f | wc -l)" -gt 0

    echo "Callable-sites calculation complete."
  >>>

  output {
    File callable_sites = glob("output/*")[0]
    File log_file = "logs/${sample_id}.callable_sites.log"
  }

  runtime {
    docker: "caoimhem/scomatic:0.3-amd64"
    cpu: 2
    memory: "8 GB"
    disks: "local-disk 100 HDD"
    preemptible: 3
    bootDiskSizeGb: 30
  }
}
