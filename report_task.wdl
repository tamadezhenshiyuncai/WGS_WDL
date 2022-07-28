version development

task report {
    input {
        File str_xlsx # from task str_hunter
        File un_rand_bam # from task un_rand_merge
        File unmap_bam  # from task unmap_merge
        File species_abundance # from task intrauterine_infection
        File cnv_calls_anno # from task exon_depth
        File sma_txt  # from task sma
        Array[File] inhouse_wgs_files # from task variant_annotation 123 jobs
        File triploid_pdf # from task triploid
        File? cnv_loh_tsv #from task cnv_loh
        File? loh_tsv #from task loh
        File anno_xls #from task batchCNV_anno
        File anno_result_xlsx # from task lumpy_anno
        File anno_xlsx # from task variant_mt
        File contamination_xls
        File nator_mops_xlsx 
        String database_path  #such as /zfsyt1/B2C_RD_P2/PMO/liufengxia/WGS/wgs_frame/wgs_clinical/bin/wgs_database/
        String go_path #such as /zfsyt1/B2C_RD_P2/USER/yansaiying/script/WGS_pipeline/software/go
        String tools_dir
        String batch_path
        String batch_number
        String project_path
        String conda_path
        Int cpu
        Float? mem
    }
    command <<<

        mkdir -p ~{batch_path}/report
        export PATH="~{tools_dir}/tools:$PATH"
        export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:~{tools_dir}/tools/lib64"
        #export GOROOT=~{go_path}/go/   #GO安装路径应该在make_json.py 脚本去解析conf文件得到
        #export GOTMPDIR=~{go_path}/tmp

        sh ~{tools_dir}/report_prepare.sh ~{project_path} ~{batch_number}
        sex=`awk '{print $NF}' ~{batch_path}/QC/sex.txt`
        if [ $sex == "F" ];then
        #    rm ~{batch_path}/annotation/part_anno/chrY_*
             cp ~{tools_dir}/etc/chrY_part1.anno ~{batch_path}/annotation/part_anno
        fi
        python3 ~{tools_dir}/creat_db.py \
            --sample ~{batch_number} \
            --path ~{batch_path} \
            --database ~{database_path}/sample.sqlite3 && \

        python3 ~{tools_dir}/wgs_var_stat.py \
            --dir ~{batch_path}/annotation/script/\
            --out ~{batch_path}/annotation/~{batch_number}.var_num.txt && \

        python3 ~{tools_dir}/wgs_excel_combine_py3_expand.py \
            --sheet0 --tie1dir ~{batch_path}/annotation/part_anno \
            --outofcat ~{batch_path}/report/~{batch_number}.wgs.tie1.xlsx \
            --suffix Tier1.spliceai.xlsx && \

        python3 ~{tools_dir}/wgs_excel_combine_py3_expand.py \
            --spliceai_cat ~{batch_path}/annotation/part_anno \
            --outfile ~{batch_path}/report/~{batch_number}.wgs.intron.xlsx \
            --suffix WGS.spliceai.xlsx && \

        ~{tools_dir}/wgs_database/pipeline_anno/anno2xlsx_v2/tier1tags/tier1tags \
            -snv ~{batch_path}/report/~{batch_number}.wgs.tie1.xlsx,~{batch_path}/report/~{batch_number}.wgs.intron_select.xlsx \
            --prefix ~{batch_path}/report/~{batch_number}.wgs.tie1_intron \
            -columns ~{tools_dir}/wgs_database/pipeline_anno/anno2xlsx_v2/template/wgs.tier1tag.title.txt && \

        ~{tools_dir}/wgs_database/pipeline_anno/anno2xlsx_v2/anno2xlsx \
            -karyotype ~{batch_path}/CNV/triploid/~{batch_number}.wgs.aneuploid.xls \
            -qc ~{batch_path}/QC/~{batch_number}.QC.txt \
            -exon ~{batch_path}/CNV/ExonDepth/all.CNV.calls.anno \
            -list ~{batch_number} \
            -prefix ~{batch_path}/report/~{batch_number}.qc.online \
            -wgs \
            -cfg ~{tools_dir}/etc/config.toml && \

        python3 \
            ~{tools_dir}/wgs_excel_combine_py3_expand.py \
            --online --tie1 ~{batch_path}/report/~{batch_number}.wgs.tie1_intron.tier1.xlsx \
            --exonCNVAndQc ~{batch_path}/report/~{batch_number}.qc.online.Tier1.xlsx \
            --outofonline ~{batch_path}/report/~{batch_number}.online.xlsx && \

        python3 \
            ~{tools_dir}/wgs_excel_combine_py3_expand.py \
            --offline --mtIntron ~{batch_path}/annotation/part_anno \
            --cnv ~{batch_path}/CNV/cnvnator/~{batch_number}.nator_mops_select.xlsx \
            --sv ~{batch_path}/CNV/lumpy/~{batch_number}.lumpy_anno.xlsx \
            --hba ~{batch_path}/CNV/batchCNV/annot.xls \
            --sma ~{batch_path}/CNV/SMA/SMA_v2.txt \
            --outfile ~{batch_path}/report/~{batch_number}.offline.xlsx \
            --mutect2 ~{batch_path}/variant/~{batch_number}.mutect2.anno.xlsx \
            --loh ~{batch_path}/CNV/loh/~{batch_number}-loh-anno.tsv \
            --str ~{batch_path}/CNV/STR/STR.~{batch_number}.xlsx \
            --triploid ~{batch_path}/CNV/triploid/~{batch_number}.xlsx \
            --Intrauterine_infection ~{batch_path}/Intrauterine_infection/coverage/~{batch_number}.sort.markdup.species_abundance.xls \
            --contamination ~{batch_path}/QC/~{batch_number}-contamination.xlsx && \

        python3 ~{tools_dir}/wgs_excel_combine_py3_expand.py \
            --qc_collect \
            --insert ~{batch_path}/QC/~{batch_number}.insert.jpg \
            --qc ~{batch_path}/QC/~{batch_number}.QC.txt \
            --var ~{batch_path}/annotation/~{batch_number}.var_num.txt \
            --qcout ~{batch_path}/report/~{batch_number}.qc_collect.xlsx && \

        python3 \
            ~{tools_dir}/link.py \
            ~{batch_path}/report/w.~{batch_number}.link.ls && \

        python3 ~{tools_dir}/wgs_excel_combine_py3_expand.py \
            --on_off_comb \
            --on_v1 ~{batch_path}/report/~{batch_number}.online.xlsx \
            --off_v1 ~{batch_path}/report/~{batch_number}.offline.xlsx && \

        source ~{conda_path}/bin/activate python3
        export PYTHONPATH=~{tools_dir}/wgs_database/pipeline_anno/exon_cnv_anno:$PYTHONPATH
        if [[ -f ~{batch_path}/CNV/triploid/~{batch_number}.xlsx ]];then
            ky_anno_python3 \
                ~{tools_dir}/CNV/wgs_cnv_anno/WGS.CnvAnno.py \
                -i ~{batch_path}/report/~{batch_number}.online.v2.xlsx \
                -g ~{batch_path}/report/~{batch_number}.offline.v2.xlsx \
                -at both \
                -st single \
                -o ~{batch_path}/report/~{batch_number}.online.v3.xlsx
        else
            cp ~{batch_path}/report/~{batch_number}.online.v2.xlsx ~{batch_path}/report/~{batch_number}.online.v3.xlsx
        fi && \
        tar --ignore-failed-read -cvhf ~{batch_path}/report/~{batch_number}.tar\
            -C ~{batch_path}/report \
            ~{batch_number}.online.v3.xlsx \
            ~{batch_number}.offline.v2.xlsx \
            ~{batch_number}.qc_collect.xlsx \
            ~{batch_number}-loh.pdf \
            ~{batch_number}-trip.pdf \
            ~{batch_number}.result_batCNV-dipin.tar.gz \
            ~{batch_number}-CNV-nator.pdf && \
        sh ~{tools_dir}/Report_TEST/demoreport.sh ~{batch_path} 
        #ls ~{batch_path}/report/*|grep -v report_~{batch_number}.sh|grep -v link.ls|grep -v rm.ls|grep -v ~{batch_number}.tar|xargs -I{}  rm -rf {} && \
        #less ~{batch_path}/report/w.~{batch_number}.rm.ls|awk '{print "rm -rf "$NF}'|sh
    >>>
    runtime {
        cpu: cpu
        memory: "~{mem}GB"
    }
}
