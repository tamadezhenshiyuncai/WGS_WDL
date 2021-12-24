version development

task variant_hc {
    input {
        String tools_dir
        String batch_path
        String chr
        String part_num # 例如 ~{part_num}
        String ref_fa
        File bqsr_bam
        File bqsr_bai
        File bam_merge_check
        String bed_split_path  #染色体拆份bed
        Int cpu
        Float? mem
    }
    command <<<
        export PATH="~{tools_dir}/tools:$PATH"
        mkdir -p ~{batch_path}/variant/script/tmp
        mkdir -p ~{batch_path}/variant/~{chr}/~{part_num}
        java -Xmx4G -XX:ParallelGCThreads=4 \
            -jar ~{tools_dir}/tools/gatk-package-4.0.11.0-local.jar HaplotypeCaller \
            --tmp-dir ~{batch_path}/variant/script/tmp \
            -ERC GVCF --correct-overlapping-quality true \
            -A BaseQuality -A MappingQuality -A QualByDepth -A MappingQualityRankSumTest -A ReadPosRankSumTest -A FisherStrand -A StrandOddsRatio -A InbreedingCoeff \
            -R ~{ref_fa} \
            -L ~{bed_split_path} \
            -I ~{bqsr_bam} \
            -O ~{batch_path}/variant/~{chr}/~{part_num}/~{chr}_~{part_num}.gvcf.gz && \
        
        java -Xmx4G -XX:ParallelGCThreads=4 \
            -jar ~{tools_dir}/tools/gatk-package-4.0.11.0-local.jar GenotypeGVCFs \
            --tmp-dir ~{batch_path}/variant/script/tmp \
            -R ~{ref_fa} \
            -V ~{batch_path}/variant/~{chr}/~{part_num}/~{chr}_~{part_num}.gvcf.gz \
            -L ~{bed_split_path} \
            -O ~{batch_path}/variant/~{chr}/~{part_num}/~{chr}_~{part_num}.vcf.gz 
    >>>
    runtime {
        cpu: cpu
        memory: "~{mem}GB"
    }
    ##此处要增加一步收集结果为scatter跑variant_annotation 作输入
    output {
        File split_gvcf = "~{batch_path}/variant/~{chr}/~{part_num}/~{chr}_~{part_num}.gvcf.gz"
        File split_gvcf_tbi = "~{batch_path}/variant/~{chr}/~{part_num}/~{chr}_~{part_num}.gvcf.gz.tbi"
        File split_vcf = "~{batch_path}/variant/~{chr}/~{part_num}/~{chr}_~{part_num}.vcf.gz"
        File split_vcf_tbi = "~{batch_path}/variant/~{chr}/~{part_num}/~{chr}_~{part_num}.vcf.gz.tbi"
        Pair[String, File] chr_vcf_pair = (chr,split_vcf)
    }
}

task concat_vcf {
    input {
        String tools_dir
        String batch_path
        String batch_number    #WDL串流程的原则是写task时尽量不要多做处理，输出是什么类型就是什么类型，所有的输入输出变体操作都在workflow块中用内建函数去处理
        Array[File] batch_part_vcfs  #上一次gather到的所有拆分后的vcf结果，通过write_lines得到的结果
        Array[File] batch_part_gvcfs #同上 gvcf的列表文件
        Array[File] batch_part_vcf_tbis
        Array[File] batch_part_gvcf_tbis
        Int cpu
        Float? mem
    }
    command <<<
        export PATH="~{tools_dir}/tools:$PATH"
        bcftools concat -a -D -q 30 -O z \
            -o ~{batch_path}/variant/~{batch_number}.wgs.vcf.gz \
            -f ~{write_lines(batch_part_vcfs)} && \
        tabix -p vcf -f ~{batch_path}/variant/~{batch_number}.wgs.vcf.gz && \
        bcftools concat -a -D -q 30 -O z \
            -o ~{batch_path}/variant/~{batch_number}.wgs.gvcf.gz \
            -f ~{write_lines(batch_part_gvcfs)} && \
        rm -rf ~{batch_path}/variant/~{batch_number}.wgs.filter.vcf.gz ~{batch_path}/variant/~{batch_number}.wgs.filter.vcf.gz.tbi && \
        mv ~{batch_path}/variant/~{batch_number}.wgs.vcf.gz ~{batch_path}/variant/~{batch_number}.wgs.filter.vcf.gz && \
        mv ~{batch_path}/variant/~{batch_number}.wgs.vcf.gz.tbi ~{batch_path}/variant/~{batch_number}.wgs.filter.vcf.gz.tbi
    >>>
    runtime {
        cpu: cpu
        memory: "~{mem}GB"
    }
    output {
        File batch_part_vcfs_list = write_lines(batch_part_vcfs)
        File batch_part_gvcfs_list = write_lines(batch_part_gvcfs)
        File concat_vcf = "~{batch_path}/variant/~{batch_number}.wgs.filter.vcf.gz"
        File concat_vcf_tbi = "~{batch_path}/variant/~{batch_number}.wgs.filter.vcf.gz.tbi"
        File concat_gvcf = "~{batch_path}/variant/~{batch_number}.wgs.gvcf.gz"
    }
}

task variant_mt {
    input {
        String tools_dir
        String batch_path
        String batch_number
        Int cpu
        Float? mem 
        File concat_vcf
        File concat_gvcf
    }
    command <<<
        export PATH="~{tools_dir}/tools:$PATH"
        python3 ~{tools_dir}/chrM_pipe.py \
            --dir ~{batch_path} \
            --outdir ~{batch_path}/variant \
            --tag ~{batch_number} \
            --config ~{tools_dir}/wgs_yantian.yaml
    >>>
    runtime {
        cpu: cpu
        memory: "~{mem}GB"
    }
    output {
        String variant_path = "~{batch_path}/variant"
        #File anno1_xlsx = "~{batch_path}/variant/~{batch_number}.mutect2.anno1.xlsx"
        File anno_xlsx = "~{batch_path}/variant/~{batch_number}.mutect2.anno.xlsx"
    }
}

task vqsr_snp {
    input {
        String tools_dir
        String batch_path
        String batch_number
        String ref_fa
        File concat_vcf
        File concat_vcf_tbi
        Int cpu
        Float? mem
    }
    command <<<
        export PATH="~{tools_dir}/tools:$PATH"
        mkdir -p ~{batch_path}/variant/script/snp_javatmp 
        java -jar ~{tools_dir}/tools/gatk-package-4.0.11.0-local.jar \
            SelectVariants \
            --tmp-dir=~{batch_path}/variant/script/snp_javatmp \
            -R ~{ref_fa} \
            --variant ~{concat_vcf} \
            -O ~{batch_path}/variant/~{batch_number}.wgs.snp.vcf \
            -select-type SNP && \
        java -jar ~{tools_dir}/tools/gatk-package-4.0.11.0-local.jar \
            VariantRecalibrator \
            --tmp-dir=~{batch_path}/variant/script/snp_javatmp \
            -R ~{ref_fa} \
            -V ~{batch_path}/variant/~{batch_number}.wgs.snp.vcf \
            --resource hapmap,known=false,training=true,truth=true,prior=15.0:~{tools_dir}/tools/bundle/hapmap_3.3.hg19.sites.vcf.gz \
            --resource omni,known=false,training=true,truth=false,prior=12.0:~{tools_dir}/tools/bundle/1000G_omni2.5.hg19.sites.vcf.gz \
            --resource 1000G,known=false,training=true,truth=false,prior=10.0:~{tools_dir}/tools/bundle/1000G_phase1.snps.high_confidence.hg19.sites.vcf.gz \
            --resource dbsnp,known=true,training=false,truth=false,prior=2.0:~{tools_dir}/tools/bundle/dbsnp_138.hg19.vcf.gz \
            -an DP -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
            -mode SNP \
            --max-gaussians 4 \
            -O ~{batch_path}/variant/~{batch_number}.wgs.snp_recal.vcf \
            --tranches-file ~{batch_path}/variant/~{batch_number}.wgs.snp_tranches.vcf && \
        java -jar ~{tools_dir}/tools/gatk-package-4.0.11.0-local.jar \
            ApplyVQSR \
            --tmp-dir=~{batch_path}/variant/script/snp_javatmp \
            -R ~{ref_fa} \
            -V ~{batch_path}/variant/~{batch_number}.wgs.snp.vcf \
            -O ~{batch_path}/variant/~{batch_number}.wgs.snp.filter.vcf \
            --truth-sensitivity-filter-level 99.0 \
            --tranches-file ~{batch_path}/variant/~{batch_number}.wgs.snp_tranches.vcf \
            --recal-file ~{batch_path}/variant/~{batch_number}.wgs.snp_recal.vcf \
            -mode SNP && \
        java -Djava.io.tmpdir=~{batch_path}/variant/script/snp_javatmp \
            -jar ~{tools_dir}/tools/gatk-package-4.0.11.0-local.jar \
            SortVcf \
            -I ~{batch_path}/variant/~{batch_number}.wgs.snp.filter.vcf \
            -O ~{batch_path}/variant/~{batch_number}.wgs.snp_sort.vcf
    >>>
    runtime {
        cpu: cpu
        memory: "~{mem}GB"
    }
    output {
        File snp_vcf = "~{batch_path}/variant/~{batch_number}.wgs.snp.vcf"
        File snp_recal_vcf = "~{batch_path}/variant/~{batch_number}.wgs.snp_recal.vcf"
        File snp_filter_vcf = "~{batch_path}/variant/~{batch_number}.wgs.snp.filter.vcf"
        File snp_sort_vcf = "~{batch_path}/variant/~{batch_number}.wgs.snp_sort.vcf"
    }
}

task vqsr_indel {
    input {
        String tools_dir
        String batch_path
        String batch_number
        String ref_fa
        File concat_vcf
        File concat_vcf_tbi
        Int cpu
        Float? mem
    }
    command <<<
        export PATH="~{tools_dir}/tools:$PATH"
        mkdir -p ~{batch_path}/variant/script/indel_javatmp
        java -jar ~{tools_dir}/tools/gatk-package-4.0.11.0-local.jar \
            SelectVariants \
            --tmp-dir=~{batch_path}/variant/script/indel_javatmp \
            -R ~{ref_fa} \
            --variant ~{concat_vcf} \
            -O ~{batch_path}/variant/~{batch_number}.wgs.indel.vcf \
            -select-type INDEL && \
        java -jar ~{tools_dir}/tools/gatk-package-4.0.11.0-local.jar \
            VariantRecalibrator \
            --tmp-dir=~{batch_path}/variant/script/indel_javatmp \
            -R ~{ref_fa} \
            -V ~{batch_path}/variant/~{batch_number}.wgs.indel.vcf \
            -resource mills,known=true,training=true,truth=true,prior=12.0:~{tools_dir}/tools/bundle/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf \
            -an DP -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
            -mode INDEL \
           --max-gaussians 4 \
           -O ~{batch_path}/variant/~{batch_number}.wgs.indel_recal.vcf \
           --tranches-file ~{batch_path}/variant/~{batch_number}.wgs.indel_tranches.vcf && \
        java -jar ~{tools_dir}/tools/gatk-package-4.0.11.0-local.jar \
            ApplyVQSR \
            --tmp-dir=~{batch_path}/variant/script/indel_javatmp \
            -R ~{ref_fa} \
            -V ~{batch_path}/variant/~{batch_number}.wgs.indel.vcf \
            -O ~{batch_path}/variant/~{batch_number}.wgs.indel.filter.vcf \
            --truth-sensitivity-filter-level 99.0 \
            --tranches-file ~{batch_path}/variant/~{batch_number}.wgs.indel_tranches.vcf \
            --recal-file ~{batch_path}/variant/~{batch_number}.wgs.indel_recal.vcf \
            -mode INDEL && \
        java -Djava.io.tmpdir=~{batch_path}/variant/script/indel_javatmp \
            -jar ~{tools_dir}/tools/gatk-package-4.0.11.0-local.jar \
            SortVcf \
            -I ~{batch_path}/variant/~{batch_number}.wgs.indel.filter.vcf \
            -O ~{batch_path}/variant/~{batch_number}.wgs.indel_sort.vcf 
    >>>
    runtime {
        cpu: cpu
        memory: "~{mem}GB"
    }
    output {
        File indel_vcf = "~{batch_path}/variant/~{batch_number}.wgs.indel.vcf"
        File indel_recal_vcf = "~{batch_path}/variant/~{batch_number}.wgs.indel_recal.vcf"
        File indel_filter_vcf = "~{batch_path}/variant/~{batch_number}.wgs.indel.filter.vcf"
        File indel_sort_vcf = "~{batch_path}/variant/~{batch_number}.wgs.indel_sort.vcf"
    }
}

task concat_snp_indel {
    input {
        String tools_dir
        String batch_path
        String batch_number
        File snp_sort_vcf
        File indel_sort_vcf
        Int cpu
        Float? mem
    }
    command <<<
        export PATH="~{tools_dir}/tools:$PATH"
        java -jar ~{tools_dir}/tools/gatk-package-4.0.11.0-local.jar \
            MergeVcfs \
            -I ~{snp_sort_vcf} \
            -I ~{indel_sort_vcf} \
            -O ~{batch_path}/variant/~{batch_number}.wgs.filter.vcf.gz && \
        tabix -p vcf -f ~{batch_path}/variant/~{batch_number}.wgs.filter.vcf.gz && \
        python3 ~{tools_dir}/wgs_vcf_add_vqsr.py \
            --vcffilter ~{batch_path}/variant/~{batch_number}.wgs.filter.vcf.gz \
            --beddir ~{tools_dir}/wgs_database/WGS_BED_SPLIT_24M_delN \
            --config ~{tools_dir}/wgs_yantian.yaml
    >>>
    runtime {
        cpu: cpu
        memory: "~{mem}GB"
    }
    output {
        File concat_filter_vcf = "~{batch_path}/variant/~{batch_number}.wgs.filter.vcf.gz"
    }
}

task concat_vcf_2 {
    input {
        String tools_dir
        String batch_path
        String batch_number
        File batch_part_vcfs_list
        File batch_part_gvcfs_list
        Int cpu
        Float? mem
    }
    command <<<
        export PATH="~{tools_dir}/tools:$PATH"
        bcftools concat -a -D -q 30 -O z \
            -o ~{batch_path}/variant/~{batch_number}.wgs.vcf.gz \
            -f ~{batch_part_vcfs_list} && \
        tabix -p vcf -f ~{batch_path}/variant/~{batch_number}.wgs.vcf.gz && \
        bcftools concat -a -D -q 30 -O z \
            -o ~{batch_path}/variant/~{batch_number}.wgs.gvcf.gz \
            -f ~{batch_part_gvcfs_list} && \
        rm -rf ~{batch_path}/variant/~{batch_number}.wgs.filter.vcf.gz ~{batch_path}/variant/~{batch_number}.wgs.filter.vcf.gz.tbi && \
        mv ~{batch_path}/variant/~{batch_number}.wgs.vcf.gz ~{batch_path}/variant/~{batch_number}.wgs.filter.vcf.gz && \
        mv ~{batch_path}/variant/~{batch_number}.wgs.vcf.gz.tbi ~{batch_path}/variant/~{batch_number}.wgs.filter.vcf.gz.tbi
    >>>
    runtime {
        cpu: cpu
        memory: "~{mem}GB"
    }
    output {
        File concat_vcf = "~{batch_path}/variant/~{batch_number}.wgs.filter.vcf.gz"
        File concat_vcf_tbi = "~{batch_path}/variant/~{batch_number}.wgs.filter.vcf.gz.tbi"
        File concat_gvcf = "~{batch_path}/variant/~{batch_number}.wgs.gvcf.gz"
    }
}

task variant_annotation {
    input {
        String tools_dir
        String batch_path
        String chr
        File sex_txt #来自 task qc_collect
        String part_num
        File split_vcf  #来自 task variant_hc
        #File qc_result_txt #来自 task qc_collect 并没在这个task中使用
        #File concat_filter_vcf #来自 task concat_snp_indel 并没在这个task中使用
        Int cpu
        Float? mem
    }
    command <<<
        unset PERLLIB
        unset PERL5LIB
        export BGICGA_HOME=~{tools_dir}/wgs_database/pipeline_anno/BGICG_Annotation
        export PATH="~{tools_dir}/tools:$PATH"
        mkdir -p ~{batch_path}/annotation/script
        mkdir -p ~{batch_path}/annotation/part_anno
        sex=`awk '{print $NF}' ~{sex_txt}` && \
        perl \
            ~{tools_dir}/wgs_database/pipeline_anno/bgi_anno/bin/bgicg_anno.pl \
            ~{tools_dir}/wgs_database/pipeline_anno/xgentic_Annotation/config_BGI59M_CG_single.2019.pl \
            ~{split_vcf} \
            -t vcf -n 9 -b 1000 -q -g $sex  2> ~{batch_path}/annotation/script/~{chr}_~{part_num}.sh.log \
            -o ~{batch_path}/annotation/part_anno/~{chr}_~{part_num}.anno && \
        
        ~{tools_dir}/wgs_database/pipeline_anno/ACMGlocalAnno/ACMGlocalAnno \
            -db ~{tools_dir}/wgs_database/pipeline_anno/ACMGlocalAnno/ACMG59.db.xlsx \
            -titleMap ~{tools_dir}/wgs_database/pipeline_anno/ACMGlocalAnno/title.txt \
            -input ~{batch_path}/annotation/part_anno/~{chr}_~{part_num}.anno \
            -output ~{batch_path}/annotation/part_anno/~{chr}_~{part_num}.anno.anno.bed.ACMG && \
        perl \
            ~{tools_dir}/wgs_database/pipeline_anno/xgentic_Annotation/update.Function.pl \
            ~{batch_path}/annotation/part_anno/~{chr}_~{part_num}.anno.anno.bed.ACMG \
            > ~{batch_path}/annotation/part_anno/~{chr}_~{part_num}.anno.anno.bed.ACMG.updateFunc && \
        
        perl \
            ~{tools_dir}/wgs_database/pipeline_anno/bgi_anno/db/hgmd/addHG19MG2bed.pl \
            ~{batch_path}/annotation/part_anno/~{chr}_~{part_num}.anno.anno.bed.ACMG.updateFunc \
            > ~{batch_path}/annotation/part_anno/~{chr}_~{part_num}.anno.anno.bed.ACMG.updateFunc.hgmd && \

        python3 \
            ~{tools_dir}/add_af.py \
            --inn ~{batch_path}/annotation/part_anno/~{chr}_~{part_num}.anno.anno.bed.ACMG.updateFunc.hgmd \
            --outf ~{batch_path}/annotation/part_anno/~{chr}_~{part_num}.anno.anno.bed.ACMG.updateFunc.hgmd.gnomAD \
            --base ~{tools_dir}/wgs_database/project_doc/gnomAD_split/~{chr}_~{part_num}.sort \
            --type gnomAD && \
        
        python3 \
            ~{tools_dir}/add_af.py \
            --inn ~{batch_path}/annotation/part_anno/~{chr}_~{part_num}.anno.anno.bed.ACMG.updateFunc.hgmd.gnomAD \
            --outf ~{batch_path}/annotation/part_anno/~{chr}_~{part_num}.anno.anno.bed.ACMG.updateFunc.hgmd.gnomAD.inhouse \
            --base ~{tools_dir}/wgs_database/project_doc/wgs_inhouse_snv/~{chr}_~{part_num}.sort \
            --type inhouse_snp && \
        
        ~{tools_dir}/wgs_database/pipeline_anno/xgentic_Annotation/anno2xlsx/anno2xlsx \
            -snv ~{batch_path}/annotation/part_anno/~{chr}_~{part_num}.anno.anno.bed.ACMG.updateFunc.hgmd.gnomAD.inhouse \
            -redis -redisAddr 10.2.1.4:6380 -wgs -acmg \
            -filter_variants ~{tools_dir}/etc/wgs.Tier1.filter_variants.txt.txt && \
        
        ~{tools_dir}/tools/py3-spliceai \
            ~{tools_dir}/spliceai_create.py \
            -b ~{batch_path}/annotation/part_anno/~{chr}_~{part_num}.anno.anno.bed.ACMG.updateFunc.hgmd.gnomAD.inhouse.WGS.xlsx \
            -p ~{batch_path}/annotation/part_anno \
            -o ~{chr}_~{part_num}.anno.anno.bed.ACMG.updateFunc.hgmd.gnomAD.inhouse.WGS \
            -c ~{tools_dir}/spliceai_yantian.yaml \
            --sheet intron --out_format excel && \

        sh ~{batch_path}/annotation/part_anno/spliceai.~{chr}_~{part_num}.anno.anno.bed.ACMG.updateFunc.hgmd.gnomAD.inhouse.WGS.sh && \
        ~{tools_dir}/tools/py3-spliceai \
            ~{tools_dir}/spliceai_merge.py \
            -b ~{batch_path}/annotation/part_anno/~{chr}_~{part_num}.anno.anno.bed.ACMG.updateFunc.hgmd.gnomAD.inhouse.WGS.xlsx \
            -p ~{batch_path}/annotation/part_anno \
            -o ~{chr}_~{part_num}.anno.anno.bed.ACMG.updateFunc.hgmd.gnomAD.inhouse.WGS \
            -c ~{tools_dir}/spliceai_yantian.yaml \
            --sheet intron --out_format excel && \
        ~{tools_dir}/tools/py3-spliceai \
            ~{tools_dir}/spliceai_create.py \
            -b ~{batch_path}/annotation/part_anno/~{chr}_~{part_num}.anno.anno.bed.ACMG.updateFunc.hgmd.gnomAD.inhouse.Tier1.xlsx \
            -p ~{batch_path}/annotation/part_anno \
            -o ~{chr}_~{part_num}.anno.anno.bed.ACMG.updateFunc.hgmd.gnomAD.inhouse.Tier1 \
            -c ~{tools_dir}/spliceai_yantian.yaml \
            --sheet filter_variants --out_format excel && \
        sh ~{batch_path}/annotation/part_anno/spliceai.~{chr}_~{part_num}.anno.anno.bed.ACMG.updateFunc.hgmd.gnomAD.inhouse.Tier1.sh && \
        ~{tools_dir}/tools/py3-spliceai \
            ~{tools_dir}/spliceai_merge.py \
            -b ~{batch_path}/annotation/part_anno/~{chr}_~{part_num}.anno.anno.bed.ACMG.updateFunc.hgmd.gnomAD.inhouse.Tier1.xlsx \
            -p ~{batch_path}/annotation/part_anno \
            -o ~{chr}_~{part_num}.anno.anno.bed.ACMG.updateFunc.hgmd.gnomAD.inhouse.Tier1 \
            -c ~{tools_dir}/spliceai_yantian.yaml \
            --sheet filter_variants --out_format excel
    >>>
    runtime {
        cpu: cpu
        memory: "~{mem}GB"
    }
    output {
        File anno = "~{batch_path}/annotation/part_anno/~{chr}_~{part_num}.anno"
        File acmg = "~{batch_path}/annotation/part_anno/~{chr}_~{part_num}.anno.anno.bed.ACMG"
        File updatefunc = "~{batch_path}/annotation/part_anno/~{chr}_~{part_num}.anno.anno.bed.ACMG.updateFunc"
        File hgmd = "~{batch_path}/annotation/part_anno/~{chr}_~{part_num}.anno.anno.bed.ACMG.updateFunc.hgmd"
        File gnomad = "~{batch_path}/annotation/part_anno/~{chr}_~{part_num}.anno.anno.bed.ACMG.updateFunc.hgmd.gnomAD"
        File ihhouse = "~{batch_path}/annotation/part_anno/~{chr}_~{part_num}.anno.anno.bed.ACMG.updateFunc.hgmd.gnomAD.inhouse"
        File inhouse_wgs = "~{batch_path}/annotation/part_anno/~{chr}_~{part_num}.anno.anno.bed.ACMG.updateFunc.hgmd.gnomAD.inhouse.WGS.xlsx"
    }
}

task exon_depth {
    input {
        String tools_dir
        String batch_path
        File qc_result_txt #来自task qc_collect
        File merge_bam
        File bam_merge_check  #来自task merge_bam_check
        String batch_number
        Int cpu
        Float? mem
    }
    command <<<
        export LD_LIBRARY_PATH=/home/wangyaoshen/local/lib:/share/app/gcc-5.2.0/lib64:/share/app/gcc-5.2.0/lib:/home/wangyaoshen/local/lib64:$LD_LIBRARY_PATH
        export PATH="~{tools_dir}/tools:$PATH"
        mkdir -p ~{batch_path}/CNV/ExonDepth
        sex=`awk '{print $NF}' ~{batch_path}/QC/sex.txt` && \
        sh ~{tools_dir}/CNV/ExonDepth/all.sh \
            ~{batch_number} \
            ~{merge_bam} \
            $sex \
            ~{batch_path}/CNV/ExonDepth
    >>>
    runtime {
        cpu: cpu
        memory: "~{mem}GB"
    }
    output {
        File bed_bam = "~{batch_path}/CNV/ExonDepth/~{batch_number}.final.merge.bam.hg19.chr.bed.bam"
        File cnv_calls_anno = "~{batch_path}/CNV/ExonDepth/all.CNV.calls.anno"
    }
}

task sma {
    input {
        String tools_dir
        Array[File] bam_list  #所有的task bqsr 任务跑完后，workflow收集结果 并生成列表文件
        Array[File] bai_list
        String batch_path
        Int cpu
        Float? mem
    }
    File bams_tmp = write_lines(bam_list)
    command <<<
        export PATH="~{tools_dir}/tools:$PATH"
        mkdir -p ~{batch_path}/CNV/SMA
        mv ~{bams_tmp} ~{bams_tmp + '.list'} 
        perl ~{tools_dir}/CNV/SMA/run_SMN_CNV_control_v2.pl \
            ~{bams_tmp + '.list'} \
            ~{tools_dir}/CNV/SMA/PP100.gene.info.bed \
            ~{tools_dir}/CNV/SMA/SMA_ex7_2_50X.control_gene.csv_tmp \
            ~{batch_path}/CNV/SMA && \
        sh ~{batch_path}/CNV/SMA/run_SMN_CNV_v2.sh
    >>>
    runtime {
        cpu: cpu
        memory: "~{mem}GB"
    }
    output {
        File sma_txt = "~{batch_path}/CNV/SMA/SMA_v2.txt"
    }
}

task cnvnator {
    input {
        String tools_dir
        Array[File] bqsr_check_files  #来自task bqsr_check  
        String cnvnator_path
        String chr
        String batch_path
        Int cpu
        Float? mem
    }
    command <<<
        export LD_LIBRARY_PATH=~{cnvnator_path}/gcc-5.2.0/lib64:$LD_LIBRARY_PATH
        export YEPPPLIBDIR=~{cnvnator_path}/yeppp-1.0.0/binaries/linux/x86_64/
        export YEPPPINCLUDEDIR=~{cnvnator_path}/yeppp-1.0.0/library/headers
        export ROOTSYS=~{cnvnator_path}/Root
        export PATH=$ROOTSYS/bin:$PATH
        export LD_LIBRARY_PATH=$ROOTSYS/lib:$LD_LIBRARY_PATH
        export MANPATH=$MANPATH:$ROOTSYS/man/man1
     
        mkdir -p ~{batch_path}/CNV/cnvnator
        sh ~{tools_dir}/CNV/CNVnator/script/cnvnator_chr.sh \
            ~{chr} \
            ~{batch_path}/CNV/cnvnator/~{chr}.root \
            ~{batch_path}/bam_chr/~{chr}.bqsr.bam \
            ~{tools_dir}/wgs_database/pipeline_anno/BGICG_Annotation_update/db/aln_db/hg19/chr \
            ~{batch_path}/CNV/cnvnator/~{chr}.CNV \
            1000 ~{tools_dir}/tools/cnvnator
    >>>
    runtime {
        cpu: cpu
        memory: "~{mem}GB"
    }
    output {
        File chr_cnv = "~{batch_path}/CNV/cnvnator/~{chr}.CNV"
    }
}

task mops {    #依赖 bqsr_check、qc_collect
    input {
        String tools_dir
        String batch_path
        String batch_number
        String chr
        Array[File] bqsr_check_files  #按染色体有多个文件，来自 task bqsr_check
        File qc_result_txt  # 按批次生成一个文件，来自 task qc_collect
        Int cpu
        Float? mem
    }
    Array[Array[String]] sex_cols = read_tsv("~{batch_path}/QC/sex.txt")
    Array[String] sex_flatten = flatten(sex_cols)
    Int idx = length(sex_flatten) -1
    String sex = sex_flatten[idx]
    command <<<
        export PATH="~{tools_dir}/tools:$PATH"
        mkdir -p ~{batch_path}/CNV/cnvnator/mops/~{chr}
        ret=`python3 ~{tools_dir}/CNV/mops/sex_link.py ~{tools_dir}/wgs_database/CNV_database/mops_control/normal19 ~{batch_path}/CNV/cnvnator/mops/~{chr} ~{batch_path}/QC/sex.txt ~{chr}`
        if [[ ~{chr} == "chrY" ]] && [[ ~{sex} == "F" ]];then
            echo -e "seqnames\tstart\tend\twidth\tstrand\tsampleName\tmedian\tmean\tCN" > "~{batch_path}/CNV/cnvnator/mops/~{chr}/~{chr}.1000.fast.cnMOPS.CNVs"
            touch "~{batch_path}/CNV/cnvnator/mops/~{chr}/~{chr}.1000.plot.pdf"
        fi              
        if [[ $ret == 1 ]];then
            exit 0
        fi
        ln -sf ~{batch_path}/bam_chr/~{chr}.bqsr.bam ~{batch_path}/CNV/cnvnator/mops/~{chr}/~{batch_number}.~{chr}.bqsr.bam
        ln -sf ~{batch_path}/bam_chr/~{chr}.bqsr.bai ~{batch_path}/CNV/cnvnator/mops/~{chr}/~{batch_number}.~{chr}.bqsr.bai
        cd ~{batch_path}/CNV/cnvnator/mops/~{chr}
        Rscript-mops \
            ~{tools_dir}/CNV/mops/cnMOPS.R \
            1000 ~{chr} fast
    >>>
    runtime {
        cpu: cpu
        memory: "~{mem}GB"
    }
    output {
        File cnmops_cnv = "~{batch_path}/CNV/cnvnator/mops/~{chr}/~{chr}.1000.fast.cnMOPS.CNVs"
        File plot_pdf = "~{batch_path}/CNV/cnvnator/mops/~{chr}/~{chr}.1000.plot.pdf"
    }
}

task cnvnator_mops_combine {
    input {
        String tools_dir
        Array[File] cnvnator_results  #来自收集task cnvnator 的不同染色体的结果
        Array[File] mops_results #来自收集task mops 的不同染色体的结果 
        String batch_path
        String batch_number
        Int cpu
        Float? mem
    }
    command <<<
        export PATH="~{tools_dir}/tools:$PATH"
        python3 ~{tools_dir}/CNV/CNVnator/script/nator_mops_annoSH.py \
            --nator ~{batch_path}/CNV/cnvnator \
            --mops ~{batch_path}/CNV/cnvnator/mops \
            --bamdir ~{batch_path} \
            --config ~{tools_dir}/wgs_yantian.yaml \
            --sample ~{batch_number} \
            --sex ~{batch_path}/QC/sex.txt \
            --outdir ~{batch_path}/CNV/cnvnator
    >>>
    runtime {
        cpu: cpu
        memory: "~{mem}GB"
    }
    output {
        File larg_cnv = "~{batch_path}/CNV/cnvnator/~{batch_number}.larg_CNV.xlsx"
        File nator_mops_xlsx = "~{batch_path}/CNV/cnvnator/~{batch_number}.nator_mops_select.xlsx"
        #File nator_plot_pdf = "~{batch_path}/CNV/cnvnator/~{batch_number}.nator_mops_select.xlsx.nator_plot.ls.data.pdf"
        #File mops_plot_pdf = "~{batch_path}/CNV/cnvnator/~{batch_number}.nator_mops_select.xlsx.mops_plot.ls.data.pdf"
    }
}
