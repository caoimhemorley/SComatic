version 1.0

workflow scomatic_mergebasecounts_v2 {
  input {
    String sample_id

    File astrocyte
    File excitatory
    File inhibitory
    File microglia
    File oligodendrocyte
    File opc
    File other

    Int nprocs = 4
  }

  call merge_basecounts {
    input:
      sample_id = sample_id,
      astrocyte = astrocyte,
      excitatory = excitatory,
      inhibitory = inhibitory,
      microglia = microglia,
      oligodendrocyte = oligodendrocyte,
      opc = opc,
      other = other,
      nprocs = nprocs
  }

  output {
    File merged_basecounts = merge_basecounts.merged_basecounts
    File log_file = merge_basecounts.log_file
  }
}

task merge_basecounts {
  input {
    String sample_id

    File astrocyte
    File excitatory
    File inhibitory
    File microglia
    File oligodendrocyte
    File opc
    File other

    Int nprocs
  }

  command <<<
    set -euo pipefail

    mkdir -p input merged logs

    # Keep inputs compressed, but rename symlinks so SComatic
    # correctly extracts the cell type from the filename.
    ln -s ~{astrocyte}       input/~{sample_id}.Astrocyte.gz
    ln -s ~{excitatory}      input/~{sample_id}.Excitatory.gz
    ln -s ~{inhibitory}      input/~{sample_id}.Inhibitory.gz
    ln -s ~{microglia}       input/~{sample_id}.Microglia.gz
    ln -s ~{oligodendrocyte} input/~{sample_id}.Oligodendrocyte.gz
    ln -s ~{opc}             input/~{sample_id}.OPC.gz
    ln -s ~{other}           input/~{sample_id}.Other.gz

    count=$(find input -name "*.gz" | wc -l)

    if [ "$count" -ne 7 ]; then
        echo "ERROR: Expected 7 input files, found $count"
        exit 1
    fi

    echo "Input files:"
    ls -lh input/

    /opt/conda/envs/SComatic/bin/python \
        /opt/SComatic/scripts/MergeCounts/MergeBaseCellCounts.py \
        --tsv_folder input \
        --outfile merged/~{sample_id}_merged_basecounts.tsv \
        > logs/~{sample_id}.MergeBaseCounts.log 2>&1

    test -s merged/~{sample_id}_merged_basecounts.tsv

    echo "Merged header:"
    grep '^#CHROM' merged/~{sample_id}_merged_basecounts.tsv || true

    pigz -p ~{nprocs} merged/~{sample_id}_merged_basecounts.tsv

    echo "Merge complete."
  >>>

  output {
    File merged_basecounts = "merged/~{sample_id}_merged_basecounts.tsv.gz"
    File log_file = "logs/~{sample_id}.MergeBaseCounts.log"
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
