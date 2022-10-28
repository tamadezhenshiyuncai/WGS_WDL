version development

import "variant_task.wdl" as variant

workflow var_anno_jobs {
    input {
        Array[Pair[String,File]] chr_vcf_pairs  #task variant_hc #flatten
        String tools_dir
        String batch_path
        File sex_txt  #task qc_collect
        File concat_filter_vcf #task concat_snp_indel
        Int anno_cpu
        Float? anno_mem
    } #此处依赖三个task,需要注意
    scatter (each_pair in chr_vcf_pairs) {
        String chr = each_pair.left
        String split_vcf = each_pair.right
        String vcf_base = basename(split_vcf,".vqsr.vcf.gz")
        String part_num = sub(vcf_base,"chr.*_","")
        call variant.variant_annotation as variant_annotation {
            input:
                tools_dir = tools_dir,
                batch_path = batch_path,
                chr = chr,
                sex_txt = sex_txt,
                part_num = part_num,
                split_vcf = split_vcf,
                cpu = anno_cpu,
                mem = anno_mem
        }
    }
    output {
        Array[File] inhouse_wgs_arr = variant_annotation.inhouse_wgs
    }
}
