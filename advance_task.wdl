version development

task triploid {
    input {
        String tools_dir
        String conda_path #此conda必须装有triploid环境，环境中装有triploid相关软件
        String batch_path
        String batch_number
        File concat_vcf  #来自task concat_vcf_2 
        File concat_vcf_tbi
        Int cpu
        Float? mem
    }
    command <<<
        source ~{conda_path}/bin/activate triploid
        mkdir -p ~{batch_path}/CNV/triploid
        ~{conda_path}/envs/triploid/bin/python ~{tools_dir}/CNV/triploid/triploid.py\
            -t 8 -d 10 -m 100 \
            ~{batch_number} \
            ~{concat_vcf} \
            ~{batch_path}/CNV/triploid/~{batch_number}
    >>>
    runtime {
        cpu: cpu
        memory: "~{mem}GB"
    }
    output {
        File triploid_xlsx = "~{batch_path}/CNV/triploid/~{batch_number}.xlsx"
        File triploid_pdf = "~{batch_path}/CNV/triploid/~{batch_number}.pdf"
    }
}

##这个地方要在workflow中给两个布尔值，看这两个loh 是否都要执行 用两个if block ；
task loh_cnv {
    input {
        File larg_cnv #来自 task cnvnator_mops_combine 该任务不一定执行
        File concat_vcf #来自 task concat_vcf_2 
        String tools_dir
        String batch_path
        String batch_number
        Int cpu
        Float? mem
    }
    command <<<
        export PATH="~{tools_dir}/tools:$PATH"
        export PYTHONPATH="~{tools_dir}/CNV/loh/:$PYTHONPATH"
        mkdir -p ~{batch_path}/CNV/loh
        python3_loh \
            -m  loh \
            --vcf ~{concat_vcf} \
            --wgs \
            --merge \
            --out ~{batch_path}/CNV/loh \
            --id ~{batch_number} \
            --cnv ~{larg_cnv}
    >>> 
    runtime {
        cpu: cpu
        memory: "~{mem}GB"
    }
    output {
        File cnv_loh_pdf = "~{batch_path}/CNV/loh/~{batch_number}-loh.pdf"
        File cnv_loh_xlsx = "~{batch_path}/CNV/loh/~{batch_number}-loh.xlsx"
        File cnv_loh_tsv = "~{batch_path}/CNV/loh/~{batch_number}-loh-anno.tsv"
    }
}

task loh {
    input {
        File concat_vcf #来自 task concat_vcf_2
        String tools_dir
        String batch_path
        String batch_number
        Int cpu
        Float? mem
    }
    command <<<
    export PATH="~{tools_dir}/tools:$PATH"
    export PYTHONPATH="~{tools_dir}/CNV/loh/:$PYTHONPATH"
    mkdir -p ~{batch_path}/CNV/loh
    python3_loh\
        -m  loh \
        --vcf ~{concat_vcf} \
        --wgs \
        --merge \
        --out ~{batch_path}/CNV/loh \
        --id ~{batch_number}
    >>>
    runtime {
        cpu: cpu
        memory: "~{mem}GB"
    }
    output {
        File loh_pdf = "~{batch_path}/CNV/loh/~{batch_number}-loh.pdf"
        File loh_xlsx = "~{batch_path}/CNV/loh/~{batch_number}-loh.xlsx"
        File loh_tsv = "~{batch_path}/CNV/loh/~{batch_number}-loh-anno.tsv"

    }
}

task str_hunter {
    input {
        String tools_dir
        String batch_path
        String batch_number
        Array[File] sort_bams #来自task align_sort 的结果gather
        Array[File] sort_bais
        String ref_fa
        File sex_txt   #来自 task collect_qc 
        File bam_merge_check #来自 task merge_bam_check 
        Int cpu
        Float? mem
    }
    #File one_bam = sort_bams[0]
    #String bam_base = basename(one_bam)
    #String one_bam_path = sub(one_bam, "~{bam_base}$","")
    #String glob_str = "~{one_bam_path + "../*/*.bam"}"
    command <<<
    export PATH="~{tools_dir}/tools:$PATH"
    mkdir -p ~{batch_path}/CNV/STR/bam
    python3 \
    ~{tools_dir}/CNV/STRhunter/together.py \
    -l 100 -e 1000 -p 4 \
    "~{batch_path + '/align/*/*/*.sort.bam'}" \
    ~{tools_dir}/CNV/STRhunter/variants.bed \
    | samtools sort -o ~{batch_path}/CNV/STR/bam/~{batch_number}.bam - && \

    samtools index ~{batch_path}/CNV/STR/bam/~{batch_number}.bam && \

    sex=`awk '{print $NF}' ~{sex_txt}` && \

    if [ $sex == "M" ];then
        sex_type="male"
    else
        sex_type="female"
    fi && \

    ExpansionHunter \
        --reference ~{ref_fa} \
        --sex $sex_type \
        --variant-catalog ~{tools_dir}/CNV/STRhunter/variants.json \
        --reads ~{batch_path}/CNV/STR/bam/~{batch_number}.bam \
        --output-prefix ~{batch_path}/CNV/STR/STR.~{batch_number}  && \

    jq \
        -r '.LocusResults[] | .Variants[] | "\(.VariantId),\(.Genotype),\(.GenotypeConfidenceInterval)"' \
        ~{batch_path}/CNV/STR/STR.~{batch_number}.json \
        |awk -F "," 'BEGIN {print "VariantId\tGenotype\tGenotypeConfidenceInterval"}{print $1"\t"$2"\t"$3}' \
        > ~{batch_path}/CNV/STR/STR.~{batch_number}.txt && \

    export PYTHONPATH=~{tools_dir}/CNV/str_anno:$PYTHONPATH

    python3 \
        -m str_anno \
        --in ~{batch_path}/CNV/STR/STR.~{batch_number}.txt \
        --out ~{batch_path}/CNV/STR/STR.~{batch_number}.tsv && \

    python3 ~{tools_dir}/csv2excel.py \
        --property Genotype,str,GenotypeConfidenceInterval,str \
        --txt ~{batch_path}/CNV/STR/STR.~{batch_number}.tsv \
        --outfile ~{batch_path}/CNV/STR/STR.~{batch_number}.xlsx \
        --add_sample
    >>>
    runtime {
        cpu: cpu
        memory: "~{mem}GB"
    }
    output {
        File str_xlsx = "~{batch_path}/CNV/STR/STR.~{batch_number}.xlsx"
    }
}

task intrauterine_infection {
    input {
        String tools_dir
        File bam_merge_check #来自 task merge_bam_check
        String batch_path
        String batch_number
        Int cpu
        Float? mem
    }
    command <<<
        export PATH="~{tools_dir}/tools:$PATH"
        mkdir -p "~{batch_path}/Intrauterine_infection"
        perl ~{tools_dir}/CNV/Intrauterine_infection/script/run_Intrauterine_infection_v2.pl \
            ~{batch_path}/bam_chr/~{batch_number}.final.merge.bam \
            ~{batch_number} \
            ~{batch_path}/Intrauterine_infection && \
        rm -rf ~{batch_path}/Intrauterine_infection/coverage
        rm -rf ~{batch_path}/Intrauterine_infection/bwa
        rm -rf ~{batch_path}/Intrauterine_infection/fq && \
        sh ~{batch_path}/Intrauterine_infection/run.sh
    >>>
    runtime {
        cpu: cpu
        memory: "~{mem}GB"
    }
    output {
        File species_abundance = "~{batch_path}/Intrauterine_infection/coverage/~{batch_number}.sort.rmdup.species_abundance.xls"
    }
}
