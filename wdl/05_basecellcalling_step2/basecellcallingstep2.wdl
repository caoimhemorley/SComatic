version 1.0

workflow scomatic_step2 {

    input {
        String sample_id
        File step1_calls

        Int min_distance = 5
        Int nprocs = 2
    }

    call basecellcalling_step2 {
        input:
            sample_id = sample_id,
            step1_calls = step1_calls,
            min_distance = min_distance,
            nprocs = nprocs
    }

    output {
        File step2_calls = basecellcalling_step2.step2_calls
        File step2_pass = basecellcalling_step2.step2_pass
        File log_file = basecellcalling_step2.log_file
    }
}

task basecellcalling_step2 {

    input {
        String sample_id
        File step1_calls

        Int min_distance
        Int nprocs
    }

    command <<<
        set -euo pipefail

        mkdir -p step2 logs resources

        echo "Running SComatic BaseCellCalling Step 2"
        echo "Sample: ~{sample_id}"

        EDITING_GZ=/opt/SComatic/RNAediting/AllEditingSites.hg38.txt.gz
        EDITING=resources/AllEditingSites.hg38.txt
        PON=/opt/SComatic/PoNs/PoN.scRNAseq.hg38.tsv

        echo "Checking SComatic resources..."

        if [ ! -s "$EDITING_GZ" ]; then
            echo "ERROR: RNA editing database not found or empty:"
            echo "  $EDITING_GZ"
            exit 1
        fi

        if [ ! -s "$PON" ]; then
            echo "ERROR: PoN not found or empty:"
            echo "  $PON"
            exit 1
        fi

        echo "RNA editing database:"
        ls -lh "$EDITING_GZ"

        echo "PoN:"
        ls -lh "$PON"

        echo "Decompressing RNA editing database..."

        gzip -dc "$EDITING_GZ" > "$EDITING"

        if [ ! -s "$EDITING" ]; then
            echo "ERROR: failed to decompress RNA editing database"
            exit 1
        fi

        echo "Decompressed RNA editing database:"
        ls -lh "$EDITING"

        echo "Checking input Step 1 file..."

        if [ ! -s "~{step1_calls}" ]; then
            echo "ERROR: Step 1 input file missing or empty:"
            echo "  ~{step1_calls}"
            exit 1
        fi

        ls -lh "~{step1_calls}"

        echo "Starting BaseCellCalling.step2.py..."

        /opt/conda/envs/SComatic/bin/python \
            /opt/SComatic/scripts/BaseCellCalling/BaseCellCalling.step2.py \
            --infile "~{step1_calls}" \
            --outfile "step2/~{sample_id}" \
            --editing "$EDITING" \
            --pon "$PON" \
            --min_distance ~{min_distance} \
            2>&1 | tee "logs/~{sample_id}.BaseCellCalling.step2.log"

        echo "Checking SComatic Step 2 output..."

        if [ ! -s "step2/~{sample_id}.calling.step2.tsv" ]; then
            echo "ERROR: SComatic did not produce expected Step 2 output:"
            echo "  step2/~{sample_id}.calling.step2.tsv"
            exit 1
        fi

        echo "Creating PASS-only output..."

        awk -F'\t' \
            '$0 ~ /^#/ || $6 == "PASS"' \
            "step2/~{sample_id}.calling.step2.tsv" \
            > "step2/~{sample_id}.calling.step2.pass.tsv"

        if [ ! -s "step2/~{sample_id}.calling.step2.pass.tsv" ]; then
            echo "ERROR: PASS-filtered Step 2 output is empty"
            exit 1
        fi

        echo "Compressing outputs..."

        pigz -p ~{nprocs} \
            "step2/~{sample_id}.calling.step2.tsv"

        pigz -p ~{nprocs} \
            "step2/~{sample_id}.calling.step2.pass.tsv"

        if [ ! -s "step2/~{sample_id}.calling.step2.tsv.gz" ]; then
            echo "ERROR: compressed Step 2 output missing"
            exit 1
        fi

        if [ ! -s "step2/~{sample_id}.calling.step2.pass.tsv.gz" ]; then
            echo "ERROR: compressed PASS output missing"
            exit 1
        fi

        echo "Final outputs:"
        ls -lh \
            "step2/~{sample_id}.calling.step2.tsv.gz" \
            "step2/~{sample_id}.calling.step2.pass.tsv.gz"

        echo "SComatic BaseCellCalling Step 2 complete."
    >>>

    output {
        File step2_calls = "step2/~{sample_id}.calling.step2.tsv.gz"
        File step2_pass = "step2/~{sample_id}.calling.step2.pass.tsv.gz"
        File log_file = "logs/~{sample_id}.BaseCellCalling.step2.log"
    }

    runtime {
        docker: "caoimhem/scomatic:0.3-amd64"
        cpu: nprocs
        memory: "16 GB"
        disks: "local-disk 300 HDD"
        preemptible: 3
        bootDiskSizeGb: 30
    }
}
