version development

import "structs.wdl" as structs
    alias Modules as Modules
    alias Resource as Resource
    alias FileFilePair as FileFilePair
    alias StringArrayFilePair as StringArrayFilePair
    alias StringArrayStringPair as StringArrayStringPair

import "align_task.wdl" as align
import "bqsr_task.wdl" as bqsr
import "lumpy_task.wdl" as lumpy
import "bed_split_workflow.wdl" as bed_split_workflow
import "bed_split24M_workflow.wdl" as bed_split24M_workflow

workflow split_chr_jobs {
    input {
        String tools_dir
        String batch_path
        String batch_number
        String ref_fa  #一定要用String 防止链接或拷贝
        Array[String] align_path_arr
        Array[StringArrayStringPair] chr_parts
        Int merge_dup_cpu
        Float? merge_dup_mem
        String omni_sites_vcf
        String high_confidence_sites_vcf
        String dbsnp_vcf
        String indels_sites_vcf
        Int bqsr_cpu
        Float? bqsr_mem
        Int check_cpu
        Float? check_mem
        Int lumpy_split_cpu
        Float? lumpy_split_mem
        Int lumpy_disco_cpu
        Float? lumpy_disco_mem
        Int part_split_cpu
        Float? part_split_mem
        Int bed_split_cpu
        Float? bed_split_mem
        Int hc_cpu
        Float? hc_mem
    }
    scatter (each_pair in chr_parts) {
        String chr = each_pair.left
        Array[String] chr_split24M_beds = each_pair.mid
        Array[String] chr_split_beds = each_pair.right
        #File raw_fix_dup_bam = "~{batch_path}/bam_chr/~{chr}.raw.dup.fix.bam"
        call align.bam_merge_dup as bam_merge_dup {
            input:
                tools_dir = tools_dir,
                align_path_arr = align_path_arr,
                chr = chr,
                batch_path = batch_path,
                cpu = merge_dup_cpu,
                mem = merge_dup_mem
        }
        call bqsr.bqsr as bqsr_t {
            input:
                tools_dir = tools_dir,
                batch_path = batch_path,
                ref_fa = ref_fa,
                raw_fix_dup_bam = bam_merge_dup.raw_fix_dup_bam,
                omni_sites_vcf = omni_sites_vcf,
                high_confidence_sites_vcf = high_confidence_sites_vcf,
                dbsnp_vcf = dbsnp_vcf,
                indels_sites_vcf = indels_sites_vcf,
                chr = chr,
                cpu = bqsr_cpu,
                mem = bqsr_mem
        }
        call bqsr.bqsr_check as bqsr_check {
            input:
                bqsr_table = bqsr_t.bqsr_table,
                bqsr_bam = bqsr_t.bqsr_bam,
                batch_path = batch_path,
                batch_number = batch_number,
                chr = chr,
                cpu = check_cpu,
                mem = check_mem
        }
        call lumpy.lumpy_chr_split as lumpy_chr_split {
            input:
                tools_dir = tools_dir,
                bqsr_check_file = bqsr_check.bqsr_check_file,
                chr = chr,
                batch_path = batch_path,
                cpu = lumpy_split_cpu,
                mem = lumpy_split_mem
        }
        call lumpy.lumpy_chr_disco as lumpy_chr_disco {
            input:
                tools_dir = tools_dir,
                bqsr_check_file = bqsr_check.bqsr_check_file,
                batch_path = batch_path,
                chr = chr,
                cpu = lumpy_disco_cpu,
                mem = lumpy_disco_mem
        }
        call bed_split_workflow.bed_split_jobs as bed_split_jobs {
            input:
                chr_split_beds = chr_split_beds,
                tools_dir = tools_dir,
                chr = chr,
                batch_path = batch_path,
                part_split_cpu = part_split_cpu,
                part_split_mem = part_split_mem,
                bam_merge_check = bqsr_check.bqsr_check_file,
                bqsr_bam = bqsr_t.bqsr_bam,
                bqsr_bai = bqsr_t.bqsr_bai,
                bed_split_cpu = bed_split_cpu,
                bed_split_mem = bed_split_mem,
        }
        call bed_split24M_workflow.bed_split24M_jobs as bed_split24M_jobs {
            input:
                chr_split24M_beds = chr_split24M_beds,
                tools_dir = tools_dir,
                chr = chr,
                batch_path = batch_path,
                bam_merge_check = bqsr_check.bqsr_check_file,
                bqsr_bam = bqsr_t.bqsr_bam,
                bqsr_bai = bqsr_t.bqsr_bai,
                ref_fa = ref_fa,
                hc_cpu = hc_cpu,
                hc_mem = hc_mem
        }
    }
    output {
        Array[File] lumpy_split_bams = lumpy_chr_split.lumpy_bam
        Array[File] lumpy_disco_bams = lumpy_chr_disco.lumpy_disco_bam
        Array[String] bam_chr_arr = bqsr_check.bam_chr
        Array[File] bam_files = bqsr_t.bqsr_bam
        Array[File] bai_files = bqsr_t.bqsr_bai
        Array[File] bqsr_check_files = bqsr_check.bqsr_check_file
        Array[Array[File]] split_vcfs = bed_split24M_jobs.split_vcfs
        Array[Array[File]] split_vcf_tbis = bed_split24M_jobs.split_vcf_tbis
        Array[Array[Pair[String,File]]] chr_vcf_pairs = bed_split24M_jobs.chr_vcf_pairs
        Array[Array[File]] split_gvcfs = bed_split24M_jobs.split_gvcfs
        Array[Array[File]] split_gvcf_tbis = bed_split24M_jobs.split_gvcf_tbis
        Array[Array[String]] split_bed_outs = bed_split_jobs.split_bed_outs
    }
}
