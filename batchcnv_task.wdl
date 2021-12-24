version development

task batchCNV {
    input {
        String tools_dir
        String batch_path
        Array[String] bam_chr  #来自task bqsr_check 最后通过write_lines写到文件中
        String batch_number
        Int cpu
        Float? mem
    }
    command <<<
        export PATH="~{tools_dir}/tools:$PATH"
        mkdir -p ~{batch_path}/CNV/batchCNV
        cd ~{batch_path}/CNV/batchCNV
        ls -d ~{batch_path}/CNV/batchCNV/*|grep -v batchCNV_~{batch_number}.sh|grep -v batchCNV_anno_~{batch_number}.sh|grep -v breakpoint|xargs -I{}   rm -rf {} && \

        perl ~{tools_dir}/CNV/batchCNV/bin/bam2sample.pl \
            ~{write_lines(bam_chr)} \
            ~{batch_path}/CNV/batchCNV/test.sample.list && \

        cp ~{tools_dir}/CNV/batchCNV/jiyinku_global.bash.cfg ~{tools_dir}/CNV/batchCNV/batCNV_control_depth/control.sample.list ~{batch_path}/CNV/batchCNV && \

        cp -fr ~{tools_dir}/CNV/batchCNV/WGS_MA/PP100_MA_bed ~{batch_path}/CNV/batchCNV/bed && \

        perl ~{tools_dir}/CNV/batchCNV/control_depth_v1/batCNV-PP100.pl \
            ~{tools_dir}/CNV/batchCNV/jiyinku_global.bash.cfg
    >>>
    runtime {
        cpu: cpu
        memory: "~{mem}GB"
    }
    output {
        File batchcnv_merged_tsv = "~{batch_path}/CNV/batchCNV/batCNV_merged_report.tsv"
    }
}

task batchCNV_breakpoint {
    input {
        String tools_dir
        String batch_path
        Array[String] bam_chr #来自task bqsr_check 最后通过write_lines写到文件中
        Int cpu
        Float? mem
    }
    command <<<
        export PATH="~{tools_dir}/tools:$PATH"
        mkdir -p ~{batch_path}/CNV/batchCNV/breakpoint
        cd ~{batch_path}/CNV/batchCNV/breakpoint
        cp ~{tools_dir}/CNV/batchCNV/WGS_MA/breakpoint-dipin-v1.1/spec.cfg ~{tools_dir}/CNV/batchCNV/WGS_MA/breakpoint-dipin-v1.1/type.txt ~{batch_path}/CNV/batchCNV/breakpoint && \
        perl ~{tools_dir}/CNV/batchCNV/bin/bam2sample.pl \
            ~{write_lines(bam_chr)} \
            ~{batch_path}/CNV/batchCNV/breakpoint/test.sample.list && \
        perl ~{tools_dir}/CNV/batchCNV/WGS_MA/breakpoint-dipin-v1.1/bin/breakpoint-dipin.pl "bash"
    >>>
    runtime {
        cpu: cpu
        memory: "~{mem}GB"
    }
    output {
        File breakpoint_tsv = "~{batch_path}/CNV/batchCNV/breakpoint/breakpoint_report.tsv"
    }
}

task batchCNV_anno {
    input {
        String tools_dir
        String batch_path
        File batchcnv_merged_tsv
        File breakpoint_tsv
        Int cpu
        Float? mem
    }
    command <<<
        export PATH="~{tools_dir}/tools:$PATH"
        cd ~{batch_path}/CNV/batchCNV && \

        perl ~{tools_dir}/CNV/batchCNV/WGS_MA/batCNV-PP100_NewBam_v2.0/bin/anno_report_v2.1.pl \
            ~{tools_dir}/CNV/batchCNV/WGS_MA/type.bed \
            ~{tools_dir}/CNV/batchCNV/WGS_MA/batCNV-PP100_NewBam_v2.0/type_overlap.bed \
            60 50 \
            ~{batch_path}/CNV/batchCNV/batCNV_merged_report.tsv \
            ~{batch_path}/CNV/batchCNV/breakpoint/breakpoint_report.tsv \
            ~{batch_path}/CNV/batchCNV/annot.xls
    >>>
    runtime {
        cpu: cpu
        memory: "~{mem}GB"
    }
    output {
        File anno_xls = "~{batch_path}/CNV/batchCNV/annot.xls"
    }
}
