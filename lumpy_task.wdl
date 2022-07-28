version development

task lumpy_chr_split {
    input {
        String tools_dir
        File bqsr_check_file  #来自task bqsr_check
        String chr
        String batch_path
        Int cpu
        Float? mem
    }
    command <<<
        export PATH="~{tools_dir}/tools:$PATH"
        mkdir -p ~{batch_path}/CNV/lumpy/tmp
        samtools \
            view -h ~{batch_path}/bam_chr/~{chr}.bqsr.bam \
            |~{tools_dir}/CNV/lumpy/scripts/extractSplitReads_BwaMem \
            -i stdin\
            |samtools view -Sb - \
            > ~{batch_path}/CNV/lumpy/tmp/~{chr}.split.bam && \
        samtools index ~{batch_path}/CNV/lumpy/tmp/~{chr}.split.bam
    >>>
    runtime {
        cpu: cpu
        memory: "~{mem}GB"
    }
    output {
        File lumpy_bam = "~{batch_path}/CNV/lumpy/tmp/~{chr}.split.bam"
    }
}

task lumpy_chr_disco {
    input {
        String tools_dir
        File bqsr_check_file #来自task bqsr_check
        String batch_path
        String chr
        Int cpu
        Float? mem
    }
    command <<<
        export PATH="~{tools_dir}/tools:$PATH"
        mkdir -p ~{batch_path}/CNV/lumpy/tmp
        samtools \
            view -b -F 1294 \
            ~{batch_path}/bam_chr/~{chr}.bqsr.bam \
            > ~{batch_path}/CNV/lumpy/tmp/~{chr}.disco.bam && \
        samtools index ~{batch_path}/CNV/lumpy/tmp/~{chr}.disco.bam
    >>>
    runtime {
        cpu: cpu
        memory: "~{mem}GB"
    }
    output {
        File lumpy_disco_bam = "~{batch_path}/CNV/lumpy/tmp/~{chr}.disco.bam"
    }
}

task final_lumpy {
    input {
        String tools_dir
        String batch_path
        String batch_number
        String conda_path  #此conda必须装有lumpy环境，环境中装有lumpy软件  
        Array[File] split_bams    #收集上一步task lumpy_chr_split 的结果 并write_lines() 生成一个文件
        Array[File] disco_bams    #收集上一步task lumpy_chr_disco 的结果 并write_lines() 生成一个文件
        File bam_merge_check  #From task merge_bam_check
        Int cpu
        Float? mem
    }
    command <<<
        export PATH="~{tools_dir}/tools:$PATH"
        samtools \
            cat -b ~{write_lines(split_bams)} \
            -o ~{batch_path}/CNV/lumpy/~{batch_number}.split.bam && \

        samtools index ~{batch_path}/CNV/lumpy/~{batch_number}.split.bam && \

        samtools \
            cat -b ~{write_lines(disco_bams)} \
            -o ~{batch_path}/CNV/lumpy/~{batch_number}.disco.bam && \

        samtools index ~{batch_path}/CNV/lumpy/~{batch_number}.disco.bam && \

        source ~{conda_path}/bin/activate lumpy    #如果用容器打包环境，这些与环境相关的语句都会在容器启动时执行 （WDL 支持 docker&singularity）
        export LD_LIBRARY_PATH=/share/app/gcc-5.2.0/lib64:$LD_LIBRARY_PATH
        export PATH=~{conda_path}/envs/lumpy/bin:$PATH
        export HEXDUMP=hexdump

        mkdir -p ~{batch_path}/CNV/lumpy/script/tmp

        ~{tools_dir}/CNV/lumpy/bin/lumpyexpress \
            -B ~{batch_path}/bam_chr/~{batch_number}.final.merge.bam \
            -S ~{batch_path}/CNV/lumpy/~{batch_number}.split.bam \
            -D ~{batch_path}/CNV/lumpy/~{batch_number}.disco.bam \
            -T ~{batch_path}/CNV/lumpy/script/tmp \
            -x ~{tools_dir}/CNV/lumpy/lumpy_exclude.bed.gz \
            -o ~{batch_path}/CNV/lumpy/~{batch_number}.lumpy.vcf
    >>>
    runtime {
        cpu: cpu
        memory: "~{mem}GB"
    }
    output {
        File final_split_bam = "~{batch_path}/CNV/lumpy/~{batch_number}.split.bam"
        File final_disco_bam = "~{batch_path}/CNV/lumpy/~{batch_number}.disco.bam"
        File final_lumpy_vcf = "~{batch_path}/CNV/lumpy/~{batch_number}.lumpy.vcf"
    }
}

task lumpy_anno {
    input {
        String tools_dir
        String batch_path
        String batch_number
        File final_lumpy_vcf #来自task final_lumpy
        Int cpu
        Float? mem
    }
    command <<<
        export PATH="~{tools_dir}/tools:$PATH"
        tmpDir="~{batch_path}/CNV/lumpy/tmp"  #这个地方不替换掉会导致下面脚本Manager()那一行报错
        export _JAVA_OPTIONS=-Djava.io.tmpdir="$tmpDir"
        export TMPDIR="$tmpDir"
        python2 ~{tools_dir}/CNV/lumpy/lumpy_annot/lumpy2gt_anno.py \
            -v ~{final_lumpy_vcf} \
            -o ~{batch_path}/CNV/lumpy/~{batch_number}.lumpy_anno.xlsx \
            -b ~{batch_path}/bam_chr/~{batch_number}.final.merge.bam \
            -d ~{batch_path}/CNV/lumpy/script/tmp \
            -n 500 -s ~{batch_number} \
            --config ~{tools_dir}/wgs.yaml
    >>>
    runtime {
        cpu: cpu
        memory: "~{mem}GB"
    }
    output {
        File anno_result_xlsx = "~{batch_path}/CNV/lumpy/~{batch_number}.lumpy_anno.xlsx"
    }
}
