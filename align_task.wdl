version development

task align {
    input {
        String tools_dir
        File clean_read_one   #此处的输入是来自上步的结果，在workflow里需要scatter 分散任务 piece开头的reads 比如 0001.fastp...
        File clean_read_two
        String batch_number
        String split_ref_fa
        String project_path
        String sample_name
        String piece 
        Int cpu
        Float? mem
    }
    command <<<
        export PATH="~{tools_dir}/tools:$PATH"
        mkdir -p ~{project_path}/~{batch_number}/align/~{sample_name}/~{piece}
        bwa \
            mem -M -Y -t 8 -T 0 \
            -R "@RG\tID:~{batch_number}\tSM:~{batch_number}\tPL:DNBSEQ\tLB:~{batch_number}" \
            ~{split_ref_fa} \
            ~{clean_read_one} \
            ~{clean_read_two} \
            | samtools view -S -b \
            -o ~{project_path}/~{batch_number}/align/~{sample_name}/~{piece}/~{sample_name}_~{piece}.align.bam -
    >>>
    runtime {
        cpu: cpu
        memory: "~{mem}GB"
    }
    output {
        File out_bam = "~{project_path}/~{batch_number}/align/~{sample_name}/~{piece}/~{sample_name}_~{piece}.align.bam"
    }
}

task align_sort {
    input {
        String tools_dir
        String batch_number
        String sample_name
        String batch_path
        String piece
        File piece_out_bam
        Int cpu
        Float? mem
    }
    command <<<
        export PATH="~{tools_dir}/tools:$PATH"
        mkdir -p ~{batch_path}/align/script/TMPDIR
        sambamba \
            sort --memory-limit 23G -l 1 -t 8 \
            --tmpdir= ~{batch_path}/align/script/TMPDIR \
            -o ~{batch_path}/align/~{sample_name}/~{piece}/~{sample_name}_~{piece}.align.sort.bam \
            ~{piece_out_bam} && \
        samtools index ~{batch_path}/align/~{sample_name}/~{piece}/~{sample_name}_~{piece}.align.sort.bam && \
        echo "align_sort done."
    >>>
    runtime {
        cpu: cpu
        memory: "~{mem}GB"
    }
    output {
        File sort_bam = "~{batch_path}/align/~{sample_name}/~{piece}/~{sample_name}_~{piece}.align.sort.bam"
        File sort_bai = "~{batch_path}/align/~{sample_name}/~{piece}/~{sample_name}_~{piece}.align.sort.bam.bai"
        #String batch_path = "~{batch_path}"
        File done_flag = stdout()
    }
}

task un_rand_merge{
    input {
        String tools_dir
        String batch_path
        Array[Array[File]] done_flag
        Int cpu
        Float? mem
    }
    command <<<
        export PATH="~{tools_dir}/tools:$PATH"
        mkdir -p ~{batch_path}/bam_chr
        python3 ~{tools_dir}/wgs_unmap_un_random_merge.py \
            --dir ~{batch_path} \
            --out ~{batch_path}/bam_chr/un_rand.bam \
            --kind un_rand --thread 5
    >>>
    runtime {
        cpu: cpu
        memory: "~{mem}GB"
    }
    output {
        File un_rand_bam = "~{batch_path}/bam_chr/un_rand.bam"
    }
}

task unmap_merge{
    input {
        String tools_dir
        String batch_path
        Array[Array[File]] done_flag
        Int cpu
        Float? mem
    }
    command <<<
        export PATH="~{tools_dir}/tools:$PATH"
        mkdir -p ~{batch_path}/bam_chr
        python3 ~{tools_dir}/wgs_unmap_un_random_merge.py \
            --dir ~{batch_path} \
            --out ~{batch_path}/bam_chr/unmap.bam \
            --kind unmap --thread 5
    >>>
    runtime {
        cpu: cpu
        memory: "~{mem}GB"
    }
    output {
        File unmap_bam = "~{batch_path}/bam_chr/unmap.bam"
    }
}

task bam_merge_dup{
    input {
        String tools_dir         #总共是三层scatter 批次号|样本名|染色体
        Array[String] align_path_arr   #此路径在scatter样本名结束后要收集成Array[String],并打印输出成文件  use write_lines() buildin function
        String chr
        String batch_path
        Int cpu
        Float? mem
    }
    command <<<
        export PATH="~{tools_dir}/tools:$PATH"
        mkdir -p ~{batch_path}/bam_chr/script/~{chr}
        python3 ~{tools_dir}/wgs_chr_bam_merge.py \
            --dirf ~{write_lines(align_path_arr)} \
            --chr ~{chr} \
            --out ~{batch_path}/bam_chr/~{chr}.raw.bam && \
        java -Xmx20G -XX:ParallelGCThreads=2 \
            -Djava.io.tmpdir=~{batch_path}/bam_chr/script/~{chr} \
            -jar ~{tools_dir}/tools/gatk-package-4.0.11.0-local.jar \
            MarkDuplicates \
            --INPUT=~{batch_path}/bam_chr/~{chr}.raw.bam \
            --OUTPUT=~{batch_path}/bam_chr/~{chr}.raw.dup.bam \
            --METRICS_FILE=~{batch_path}/bam_chr/~{chr}.raw.dup.bam.metrics \
            --ASSUME_SORTED --VALIDATION_STRINGENCY SILENT --CREATE_INDEX true \
            --DUPLICATE_SCORING_STRATEGY TOTAL_MAPPED_REFERENCE_LENGTH && \
        java -Xmx20G  -XX:ParallelGCThreads=2 \
            -Djava.io.tmpdir=~{batch_path}/bam_chr/script/~{chr} \
            -jar ~{tools_dir}/tools/gatk-package-4.0.11.0-local.jar \
            FixMateInformation --VALIDATION_STRINGENCY SILENT \
            -I ~{batch_path}/bam_chr/~{chr}.raw.dup.bam \
            -O ~{batch_path}/bam_chr/~{chr}.raw.dup.fix.bam && \
        samtools index -@ 2 ~{batch_path}/bam_chr/~{chr}.raw.dup.fix.bam
    >>>
    runtime {
        cpu: cpu
        memory: "~{mem}GB"
    }
    output {
        File raw_bam = "~{batch_path}/bam_chr/~{chr}.raw.bam"
        File raw_dup_bam = "~{batch_path}/bam_chr/~{chr}.raw.dup.bam"
        File raw_fix_dup_bam = "~{batch_path}/bam_chr/~{chr}.raw.dup.fix.bam"
    }
}
