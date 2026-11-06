version 1.0

workflow Step2PassUMI {
  input {
    String sample_id
    File step2_calls

    File astrocyte_bam
    File astrocyte_bai
    File excitatory_bam
    File excitatory_bai
    File inhibitory_bam
    File inhibitory_bai
    File microglia_bam
    File microglia_bai
    File oligodendrocyte_bam
    File oligodendrocyte_bai
    File opc_bam
    File opc_bai
    File other_bam
    File other_bai

    Int max_records = 1000000000
    Int min_bq = 20
    Int min_mq = 255
    Int max_depth = 8000
  }

  call annotate_library_umi {
    input:
      sample_id = sample_id,
      step2_calls = step2_calls,
      astrocyte_bam = astrocyte_bam,
      astrocyte_bai = astrocyte_bai,
      excitatory_bam = excitatory_bam,
      excitatory_bai = excitatory_bai,
      inhibitory_bam = inhibitory_bam,
      inhibitory_bai = inhibitory_bai,
      microglia_bam = microglia_bam,
      microglia_bai = microglia_bai,
      oligodendrocyte_bam = oligodendrocyte_bam,
      oligodendrocyte_bai = oligodendrocyte_bai,
      opc_bam = opc_bam,
      opc_bai = opc_bai,
      other_bam = other_bam,
      other_bai = other_bai,
      max_records = max_records,
      min_bq = min_bq,
      min_mq = min_mq,
      max_depth = max_depth
  }

  output {
    File celltype_summary = annotate_library_umi.celltype_summary
    File variant_summary = annotate_library_umi.variant_summary
    File molecule_evidence = annotate_library_umi.molecule_evidence
    File log_file = annotate_library_umi.log_file
  }
}

task annotate_library_umi {
  input {
    String sample_id
    File step2_calls

    File astrocyte_bam
    File astrocyte_bai
    File excitatory_bam
    File excitatory_bai
    File inhibitory_bam
    File inhibitory_bai
    File microglia_bam
    File microglia_bai
    File oligodendrocyte_bam
    File oligodendrocyte_bai
    File opc_bam
    File opc_bai
    File other_bam
    File other_bai

    Int max_records
    Int min_bq
    Int min_mq
    Int max_depth
  }

  command <<<
    set -euo pipefail

    mkdir -p bam_inputs umi logs

    ln -s "~{astrocyte_bam}"       "bam_inputs/~{sample_id}.Astrocyte.bam"
    ln -s "~{astrocyte_bai}"       "bam_inputs/~{sample_id}.Astrocyte.bam.bai"

    ln -s "~{excitatory_bam}"      "bam_inputs/~{sample_id}.Excitatory.bam"
    ln -s "~{excitatory_bai}"      "bam_inputs/~{sample_id}.Excitatory.bam.bai"

    ln -s "~{inhibitory_bam}"      "bam_inputs/~{sample_id}.Inhibitory.bam"
    ln -s "~{inhibitory_bai}"      "bam_inputs/~{sample_id}.Inhibitory.bam.bai"

    ln -s "~{microglia_bam}"       "bam_inputs/~{sample_id}.Microglia.bam"
    ln -s "~{microglia_bai}"       "bam_inputs/~{sample_id}.Microglia.bam.bai"

    ln -s "~{oligodendrocyte_bam}" "bam_inputs/~{sample_id}.Oligodendrocyte.bam"
    ln -s "~{oligodendrocyte_bai}" "bam_inputs/~{sample_id}.Oligodendrocyte.bam.bai"

    ln -s "~{opc_bam}"             "bam_inputs/~{sample_id}.OPC.bam"
    ln -s "~{opc_bai}"             "bam_inputs/~{sample_id}.OPC.bam.bai"

    ln -s "~{other_bam}"           "bam_inputs/~{sample_id}.Other.bam"
    ln -s "~{other_bai}"           "bam_inputs/~{sample_id}.Other.bam.bai"

    python -u /opt/SComatic/scripts/step2_pass_umi_library.py \
      --step2 "~{step2_calls}" \
      --bam-dir bam_inputs \
      --sample "~{sample_id}" \
      --out-prefix "umi/~{sample_id}" \
      --max-records ~{max_records} \
      --min-bq ~{min_bq} \
      --min-mq ~{min_mq} \
      --max-depth ~{max_depth} \
      > "logs/~{sample_id}.umi_library.log" 2>&1

    test -s "umi/~{sample_id}.umi_celltype_summary.tsv"
    test -s "umi/~{sample_id}.umi_variant_summary.tsv"
    test -s "umi/~{sample_id}.umi_molecules.tsv"

    echo "Library-wide UMI annotation complete."
  >>>

  output {
    File celltype_summary = "umi/${sample_id}.umi_celltype_summary.tsv"
    File variant_summary = "umi/${sample_id}.umi_variant_summary.tsv"
    File molecule_evidence = "umi/${sample_id}.umi_molecules.tsv"
    File log_file = "logs/${sample_id}.umi_library.log"
  }

  runtime {
    docker: "caoimhem/scomatic:0.3-amd64"
    cpu: 1
    memory: "16 GB"
    disks: "local-disk 100 HDD"
    preemptible: 3
    bootDiskSizeGb: 30
  }
}
