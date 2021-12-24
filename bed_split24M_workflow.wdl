version development

import "variant_task.wdl" as variant

workflow bed_split24M_jobs {
    input {
        Array[String] chr_split24M_beds
        String tools_dir
        String chr
        String batch_path
        File bam_merge_check
        File bqsr_bam
        File bqsr_bai
        String ref_fa
        Int hc_cpu
        Float? hc_mem 
    }
    scatter (each_bed in chr_split24M_beds) {
        String bed_base = basename(each_bed,".bed")
        String part_num = sub(bed_base,"chr.*\\.","")
        call variant.variant_hc as variant_hc {
            input:
                tools_dir = tools_dir,
                batch_path = batch_path,
                chr = chr,
                part_num = part_num,
                ref_fa = ref_fa,
                bqsr_bam = bqsr_bam,
                bqsr_bai = bqsr_bai,
                bam_merge_check = bam_merge_check,
                bed_split_path = each_bed,
                cpu = hc_cpu,
                mem = hc_mem
        }
    }
    output {
        Array[File] split_vcfs = variant_hc.split_vcf
        Array[File] split_vcf_tbis = variant_hc.split_vcf_tbi
        Array[File] split_gvcfs = variant_hc.split_gvcf
        Array[File] split_gvcf_tbis = variant_hc.split_gvcf_tbi
        Array[Pair[String,File]] chr_vcf_pairs = variant_hc.chr_vcf_pair
    }
}
