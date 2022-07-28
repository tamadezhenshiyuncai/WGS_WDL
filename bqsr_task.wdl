version development

task bqsr{
    input {
        String tools_dir
        String batch_path
        String ref_fa
        File raw_fix_dup_bam
        String omni_sites_vcf
        String high_confidence_sites_vcf
        String dbsnp_vcf
        String indels_sites_vcf
        String chr
        Int cpu
        Float? mem
    }
    command <<<
        export PATH="~{tools_dir}/tools:$PATH"
        mkdir -p ~{batch_path}/bam_chr/script/~{chr}_tmp
        java -jar ~{tools_dir}/tools/gatk-package-4.0.11.0-local.jar \
        BaseRecalibrator \
        -R ~{ref_fa} \
        -I ~{raw_fix_dup_bam} \
        --tmp-dir ~{batch_path}/bam_chr/script/~{chr}_tmp \
        --known-sites ~{omni_sites_vcf} \
        --known-sites ~{high_confidence_sites_vcf} \
        --known-sites ~{dbsnp_vcf} \
        --known-sites ~{indels_sites_vcf} \
        -O ~{batch_path}/bam_chr/script/~{chr}.bqsr.table && \
        java -Djava.io.tmpdir=~{batch_path}/bam_chr/script/~{chr}_tmp \
        -jar ~{tools_dir}/tools/gatk-package-4.0.11.0-local.jar \
        ApplyBQSR \
        -I  ~{raw_fix_dup_bam} \
        -bqsr ~{batch_path}/bam_chr/script/~{chr}.bqsr.table \
        -O ~{batch_path}/bam_chr/~{chr}.bqsr.bam
    >>>
    runtime {
        cpu: cpu
        memory: "~{mem}GB"
        maxRetries: 3
    }
    output {
        File bqsr_table = "~{batch_path}/bam_chr/script/~{chr}.bqsr.table"
        File bqsr_bam = "~{batch_path}/bam_chr/~{chr}.bqsr.bam"
        File bqsr_bai = "~{batch_path}/bam_chr/~{chr}.bqsr.bai"
    }
}

task bqsr_check{
    input {
        File bqsr_table
        File bqsr_bam
        String batch_path
        String batch_number
        String chr
        Int cpu
        Float? mem
    }
    command <<<
        echo "The analysis of bqsr had done! This is a startpoint!" && \
        touch ~{batch_path}/bam_chr/script/~{batch_number}.~{chr}.bqsr_done.check && \
        sleep 10s
    >>>
    runtime {
        cpu: cpu
        memory: "~{mem}GB"
    }
    output {
        File bqsr_check_file = "~{batch_path}/bam_chr/script/~{batch_number}.~{chr}.bqsr_done.check"
        String bam_chr = "~{batch_path}/bam_chr"
    }
}

task all_chr_merge_dup{
    input {
        String tools_dir
        Array[File] bqsr_check_files
        Array[File] bam_list         #这一文件建议直接由make_json.py脚本根据输入配置成后，在workflow中用write_lines()函数生成配置文件
        String batch_number
        String batch_path
        Int cpu
        Float? mem
    }
    command <<<
        export PATH="~{tools_dir}/tools:$PATH"
        samtools cat \
            -b ~{write_lines(bam_list)} \
            -o ~{batch_path}/bam_chr/~{batch_number}.final.merge.bam && \
        samtools index -@ 10 ~{batch_path}/bam_chr/~{batch_number}.final.merge.bam
    >>>
    runtime {
        cpu: cpu
        memory: "~{mem}GB"
        maxRetries: 3
    }
    output {
        File merge_bam = "~{batch_path}/bam_chr/~{batch_number}.final.merge.bam"
    }
}

task merge_bam_check {
    input {
        String batch_path
        String batch_number
        File merge_bam
        Int cpu
        Float? mem
    }
    command <<<
        echo "The analysis of all_chr_merge_dup had done! This is a startpoint!" && \
        touch ~{batch_path}/bam_chr/~{batch_number}.merge_done.check && \
        sleep 10s
    >>>
    runtime {
        cpu: cpu
        memory: "~{mem}GB"
    }
    output {
        File bam_merge_check = "~{batch_path}/bam_chr/~{batch_number}.merge_done.check"
    }
}

task bam_chr_part_split {
    input {
        String tools_dir
        String batch_path
        String bed_split_path
        File bam_merge_check
        File bqsr_bam # bqsq那一步生成的bam 因此这与那个任务应该是要在一个scatter中的   
        File bqsr_bai
        String chr   #这几个数据都要在workflow中提供是要处理好传进来的，要让task部分的输入输出尽量简单
        String part_num #比如part1
        Int cpu
        Float? mem
    }
    command <<<
        export PATH="~{tools_dir}/tools:$PATH"
        mkdir -p ~{batch_path}/QC/bam_chr_part_split/~{chr}/~{part_num}
        sambamba view -f bam -h \
        -o ~{batch_path}/QC/bam_chr_part_split/~{chr}/~{part_num}/~{chr}_~{part_num}.bam \
        -L ~{bed_split_path} \
        ~{bqsr_bam} && \
        sambamba index ~{batch_path}/QC/bam_chr_part_split/~{chr}/~{part_num}/~{chr}_~{part_num}.bam
    >>>
    runtime {
        cpu: cpu
        memory: "~{mem}GB"
    }
    output {
        File split_bam = "~{batch_path}/QC/bam_chr_part_split/~{chr}/~{part_num}/~{chr}_~{part_num}.bam"
    }
}

task qc_bed_split {
    input {
        String tools_dir
        String batch_path
        String bed_split_path
        String chr
        String part_num
        File split_bam
        Int cpu
        Float? mem
    }
    command <<<
        export PATH="~{tools_dir}/tools:$PATH"
        mkdir -p ~{batch_path}/QC/QC_BED_SPLIT/~{chr}/~{part_num}
        bamdst -p ~{bed_split_path} -o ~{batch_path}/QC/QC_BED_SPLIT/~{chr}/~{part_num} \
        ~{split_bam} || true
    >>>
    runtime {
        cpu: cpu
        memory: "~{mem}GB"
    }
    output {
        String split_bed = "~{batch_path}/QC/QC_BED_SPLIT/~{chr}/~{part_num}"
        String filter_path = "~{batch_path}/filter"  
        String qc_path = "~{batch_path}/QC"     
    }
}

task qc_collect {
    input {
        String batch_path
        String batch_number
        String tools_dir
        String sex
        Array[Array[String]] split_bed_outs
        Int cpu
        Float? mem
    }
    command <<<
        export PATH="~{tools_dir}/tools:$PATH"
        mkdir -p ~{batch_path}/filter
        mkdir -p ~{batch_path}/QC
        python2 ~{tools_dir}/wgs_qc_collect.py \
            --filter ~{batch_path}/filter \
            --bamqc ~{batch_path}/QC \
            --sex ~{sex} \
            --config ~{tools_dir}/wgs_yantian.yaml \
            > ~{batch_path}/QC/~{batch_number}.QC.txt 
    >>>
    runtime {
        cpu: cpu
        memory: "~{mem}GB"
    }
    output {
        File qc_result_txt = "~{batch_path}/QC/~{batch_number}.QC.txt"
        File sex_txt = "~{batch_path}/QC/sex.txt"
    }
}
