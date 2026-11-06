version 1.0

workflow SplitBamWorkflow {

    input {
        File bam
        File bai
        File metadata_tsv
        String sample_id

        Int min_MQ = 255
        Int max_nM = 5
        Int max_NH = 1
        Int n_trim = 0

        Int cpu = 4
        Int memory_gb = 16
    }

    call SplitBam {
        input:
            bam = bam,
            bai = bai,
            metadata_tsv = metadata_tsv,
            sample_id = sample_id,
            min_MQ = min_MQ,
            max_nM = max_nM,
            max_NH = max_NH,
            n_trim = n_trim,
            cpu = cpu,
            memory_gb = memory_gb
    }

    output {
        File astrocyte_bam = SplitBam.astrocyte_bam
        File astrocyte_bai = SplitBam.astrocyte_bai

        File excitatory_bam = SplitBam.excitatory_bam
        File excitatory_bai = SplitBam.excitatory_bai

        File inhibitory_bam = SplitBam.inhibitory_bam
        File inhibitory_bai = SplitBam.inhibitory_bai

        File microglia_bam = SplitBam.microglia_bam
        File microglia_bai = SplitBam.microglia_bai

        File oligodendrocyte_bam = SplitBam.oligodendrocyte_bam
        File oligodendrocyte_bai = SplitBam.oligodendrocyte_bai

        File opc_bam = SplitBam.opc_bam
        File opc_bai = SplitBam.opc_bai

        File other_bam = SplitBam.other_bam
        File other_bai = SplitBam.other_bai

        File report = SplitBam.report
    }
}

task SplitBam {

    input {
        File bam
        File bai
        File metadata_tsv
        String sample_id

        Int min_MQ
        Int max_nM
        Int max_NH
        Int n_trim

        Int cpu
        Int memory_gb
    }

    command <<<
        set -euo pipefail

        mkdir -p output

        python3 /opt/SComatic/scripts/SplitBam/SplitBamCellTypes.py \
            --bam ~{bam} \
            --meta ~{metadata_tsv} \
            --id ~{sample_id} \
            --min_MQ ~{min_MQ} \
            --max_nM ~{max_nM} \
            --max_NH ~{max_NH} \
            --n_trim ~{n_trim} \
            --outdir output
    >>>

    output {
        File astrocyte_bam = "output/~{sample_id}.Astrocyte.bam"
        File astrocyte_bai = "output/~{sample_id}.Astrocyte.bam.bai"

        File excitatory_bam = "output/~{sample_id}.Excitatory.bam"
        File excitatory_bai = "output/~{sample_id}.Excitatory.bam.bai"

        File inhibitory_bam = "output/~{sample_id}.Inhibitory.bam"
        File inhibitory_bai = "output/~{sample_id}.Inhibitory.bam.bai"

        File microglia_bam = "output/~{sample_id}.Microglia.bam"
        File microglia_bai = "output/~{sample_id}.Microglia.bam.bai"

        File oligodendrocyte_bam = "output/~{sample_id}.Oligodendrocyte.bam"
        File oligodendrocyte_bai = "output/~{sample_id}.Oligodendrocyte.bam.bai"

        File opc_bam = "output/~{sample_id}.OPC.bam"
        File opc_bai = "output/~{sample_id}.OPC.bam.bai"

        File other_bam = "output/~{sample_id}.Other.bam"
        File other_bai = "output/~{sample_id}.Other.bam.bai"

        File report = "output/~{sample_id}.report.txt"
    }

    runtime {
        docker: "caoimhem/scomatic:0.3-amd64"
        cpu: cpu
        memory: "16 GB"
        disks: "local-disk 250 HDD"
        preemptible: 3
        bootDiskSizeGb: 30
    }
}
