version development

import "bqsr_task.wdl" as bqsr
#import "variant_task.wdl" as variant

workflow bed_split_jobs {
    input {
        Array[String] chr_split_beds
        String tools_dir
        String  chr
        String batch_path
        Int part_split_cpu
        Float? part_split_mem
        File bam_merge_check
        File bqsr_bam
        File bqsr_bai
        Int bed_split_cpu
        Float? bed_split_mem
        #String ref_fa
    }
    scatter (each_bed in chr_split_beds) {
        String bed_base = basename(each_bed,".bed")
        String part_num = sub(bed_base,"chr.*\\.","")
        call bqsr.bam_chr_part_split as bam_chr_part_split {
            input:
                tools_dir = tools_dir,
                batch_path = batch_path,
                bed_split_path = each_bed,
                bam_merge_check = bam_merge_check,
                bqsr_bam = bqsr_bam,
                bqsr_bai = bqsr_bai,
                chr = chr,
                part_num = part_num,
                cpu = part_split_cpu,
                mem = part_split_mem
        }
        call bqsr.qc_bed_split as qc_bed_split {
            input:
                tools_dir = tools_dir,
                batch_path = batch_path,
                bed_split_path = each_bed,
                chr = chr,
                part_num = part_num,
                split_bam = bqsr_bam,
                cpu = bed_split_cpu,
                mem = bed_split_mem
        }
    }
    output {
        #Array[File] split_vcfs = variant_hc.split_vcf
        #Array[File] split_gvcfs = variant_hc.split_gvcf
        #Array[File] inhouse_wgs_files = variant_annotation.inhouse_wgs
        Array[String] split_bed_outs = qc_bed_split.split_bed
    }
}
