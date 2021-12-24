version development

import "structs.wdl" as structs
    alias Modules as Modules
    alias Resource as Resource
    alias FileFilePair as FileFilePair
    alias StringArrayFilePair as StringArrayFilePair
    alias StringArrayStringPair as StringArrayStringPair
#import "qc_task.wdl" as qc
import "align_task.wdl" as align
import "bqsr_task.wdl" as bqsr
import "variant_task.wdl" as variant
import "batchcnv_task.wdl" as batchcnv
import "lumpy_task.wdl" as lumpy
import "advance_task.wdl" as advance
import "report_task.wdl" as report
import "paired_reads_workflow.wdl" as paired_reads_workflow
import "chr_split_workflow.wdl" as chr_split_workflow
import "variant_anno_workflow.wdl" as variant_anno_workflow

workflow wgs_pip {
    input {
        Modules resources
        String tools_dir
        #Array[String] batch_nums
        String project_path
        String ref_fa
        Array[StringArrayFilePair] reads_info
        Array[String] chrs
        Array[String] noM_chrs
        String omni_sites_vcf
        String high_confidence_sites_vcf
        String dbsnp_vcf
        String indels_sites_vcf
        Array[StringArrayStringPair] chr_parts
        String sex
        String cnvnator_path
        String conda_lumpy_env
        String conda_triploid_env
        Boolean cnvloh_isrun
        String database_path
        String go_path
    }

    meta {
        author: "fuxiangke"
        email: "fuxiangke@genomics.cn"
    }

    Resource resource_fastp = resources.filter_split
    Resource resource_bqsr = resources.bqsr
    Resource resource_bqsr_check = resources.bqsr_check
    Resource resource_all_check_file = resources.all_chr_merge_dup
    Resource resource_merge_bam_check = resources.merge_bam_check
    Resource resource_bam_chr_part_split = resources.bam_chr_part_split
    Resource resource_qc_bed_split = resources.qc_bed_split
    Resource resource_qc_collect = resources.qc_collect
    Resource resource_variant_hc = resources.variant_hc
    Resource resource_align = resources.align
    Resource resource_align_sort = resources.align_sort
    Resource resource_un_rand_merge = resources.un_rand_merge
    Resource resource_unmap_merge = resources.unmap_merge
    Resource resource_bam_merge_dup = resources.bam_merge_dup
    Resource resource_concat_vcf = resources.concat_vcf
    Resource resource_variant_mt = resources.variant_mt
    Resource resource_vqsr_snp = resources.vqsr_snp
    Resource resource_vqsr_indel = resources.vqsr_indel
    Resource resource_concat_snp_indel = resources.concat_snp_indel
    Resource resource_concat_vcf_2 = resources.concat_vcf_2
    Resource resource_variant_annotation = resources.variant_annotation
    Resource resource_exon_depth = resources.exon_depth
    Resource resource_sma = resources.sma
    Resource resource_cnvnator = resources.cnvnator
    Resource resource_mops = resources.mops
    Resource resource_cnvnator_mops_combine = resources.cnvnator_mops_combine
    Resource resource_lumpy_chr_split = resources.lumpy_chr_split
    Resource resource_lumpy_chr_disco = resources.lumpy_chr_disco
    Resource resource_final_lumpy = resources.final_lumpy
    Resource resource_lumpy_anno = resources.lumpy_anno
    Resource resource_batchCNV = resources.batchCNV
    Resource resource_batchCNV_breakpoint = resources.batchCNV_breakpoint
    Resource resource_batchCNV_anno = resources.batchCNV_anno
    Resource resource_triploid = resources.triploid
    Resource resource_loh = resources.loh
    Resource resource_str_hunter = resources.str_hunter
    Resource resource_intrauterine_infection = resources.intrauterine_infection
    Resource resource_report = resources.report
        
    scatter (each_batch_info in reads_info){
        String batch_number = each_batch_info.left
        String batch_path = "~{project_path}/~{batch_number}"
        call paired_reads_workflow.paired_reads_fastp as qc_align {
            input:
                batch_path = batch_path,
                batch_number = batch_number,
                tools_dir = tools_dir,
                project_path = project_path,
                ref_fa = ref_fa,
                qc_cpu = resource_fastp.cpu,
                qc_mem  = resource_fastp.mem,
                align_cpu = resource_align.cpu,
                align_mem = resource_align.mem,
                sort_cpu = resource_align_sort.cpu,
                sort_mem = resource_align_sort.mem,
                batch_paired_reads = each_batch_info.right
        }
        #######
        Array[String] align_path_arr = qc_align.align_path_arr
        Array[Array[String]] done_flag = qc_align.done_flag
        call align.un_rand_merge as un_rand_merge {
            input:
                tools_dir = tools_dir,
                batch_path = batch_path,
                done_flag = done_flag,
                cpu = resource_un_rand_merge.cpu,
                mem = resource_un_rand_merge.mem
        }
        call align.unmap_merge as unmap_merge {
            input:
                tools_dir = tools_dir,
                batch_path = batch_path,
                done_flag = done_flag,
                cpu = resource_unmap_merge.cpu,
                mem = resource_unmap_merge.mem
        }
        call chr_split_workflow.split_chr_jobs as split_chr_jobs {
            input:
                tools_dir = tools_dir,
                batch_path = batch_path,
                batch_number = batch_number,
                ref_fa = ref_fa,
                align_path_arr = align_path_arr,
                chr_parts = chr_parts,
                merge_dup_cpu = resource_bam_merge_dup.cpu,
                merge_dup_mem = resource_bam_merge_dup.mem,
                omni_sites_vcf = omni_sites_vcf,
                high_confidence_sites_vcf = high_confidence_sites_vcf,
                dbsnp_vcf = dbsnp_vcf,
                indels_sites_vcf = indels_sites_vcf,
                bqsr_cpu = resource_bqsr.cpu,
                bqsr_mem = resource_bqsr.mem,
                check_cpu = resource_bqsr_check.cpu,
                check_mem = resource_bqsr_check.mem,
                lumpy_split_cpu = resource_lumpy_chr_split.cpu,
                lumpy_split_mem = resource_lumpy_chr_split.mem,
                lumpy_disco_cpu = resource_lumpy_chr_disco.cpu,
                lumpy_disco_mem = resource_lumpy_chr_disco.mem,
                part_split_cpu = resource_bam_chr_part_split.cpu,
                part_split_mem = resource_bam_chr_part_split.mem,
                bed_split_cpu = resource_qc_bed_split.cpu,
                bed_split_mem = resource_qc_bed_split.mem,
                hc_cpu = resource_variant_hc.cpu,
                hc_mem = resource_variant_hc.mem
        }
        ##跑vcf注释部分
        ##flatten###这一步的结果只有task report调用到
        Array[Pair[String,File]] vcf_pairs_arr = flatten(split_chr_jobs.chr_vcf_pairs)
        call variant_anno_workflow.var_anno_jobs as var_anno_jobs {
            input:
                chr_vcf_pairs = vcf_pairs_arr,
                tools_dir = tools_dir,
                batch_path = batch_path,
                sex_txt = qc_collect.sex_txt,
                concat_filter_vcf = concat_snp_indel.concat_filter_vcf,
                anno_cpu = resource_variant_annotation.cpu,
                anno_mem = resource_variant_annotation.mem
        }
        ##flatten
        Array[Array[File]] split_vcfs = split_chr_jobs.split_vcfs
        Array[Array[File]] split_vcf_tbis = split_chr_jobs.split_vcf_tbis
        Array[Array[File]] split_gvcfs = split_chr_jobs.split_gvcfs
        Array[Array[File]] split_gvcf_tbis = split_chr_jobs.split_gvcf_tbis
        Array[File] split_vcfs_flatten = flatten(split_vcfs)
        Array[File] split_vcf_tbis_flatten = flatten(split_vcf_tbis)
        Array[File] split_gvcfs_flatten = flatten(split_gvcfs)
        Array[File] split_gvcf_tbis_flatten = flatten(split_gvcf_tbis)
        Array[File] lumpy_split_bams = split_chr_jobs.lumpy_split_bams
        Array[File] lumpy_disco_bams = split_chr_jobs.lumpy_disco_bams
        Array[String] bam_chr_arr = split_chr_jobs.bam_chr_arr
        Array[String] one_bam_chr = [bam_chr_arr[0]]
        Array[File] inhouse_wgs_files = var_anno_jobs.inhouse_wgs_arr
        ##batchCNV
        call batchcnv.batchCNV as batchCNV {
            input:
                tools_dir = tools_dir,
                batch_path = batch_path,
                bam_chr = one_bam_chr,
                batch_number = batch_number,
                cpu = resource_batchCNV.cpu,
                mem = resource_batchCNV.mem
        }
        call batchcnv.batchCNV_breakpoint as batchCNV_breakpoint {
            input:
                tools_dir = tools_dir,
                batch_path = batch_path,
                bam_chr = one_bam_chr,
                cpu = resource_batchCNV_breakpoint.cpu,
                mem = resource_batchCNV_breakpoint.mem
        }
        call batchcnv.batchCNV_anno as batchCNV_anno {
            input:
                tools_dir = tools_dir,
                batch_path = batch_path,
                batchcnv_merged_tsv = batchCNV.batchcnv_merged_tsv,
                breakpoint_tsv = batchCNV_breakpoint.breakpoint_tsv,
                cpu = resource_batchCNV_anno.cpu,
                mem = resource_batchCNV_anno.mem
        }
        call lumpy.final_lumpy as final_lumpy {
            input:
                tools_dir = tools_dir,
                batch_path = batch_path,
                batch_number = batch_number,
                conda_path = conda_lumpy_env,
                split_bams = lumpy_split_bams,
                disco_bams = lumpy_disco_bams,
                bam_merge_check = merge_bam_check.bam_merge_check,
                cpu = resource_final_lumpy.cpu,
                mem = resource_final_lumpy.mem
        }
        call lumpy.lumpy_anno as lumpy_anno {
            input:
                tools_dir = tools_dir,
                batch_path = batch_path,
                batch_number = batch_number,
                final_lumpy_vcf = final_lumpy.final_lumpy_vcf,
                cpu = resource_lumpy_anno.cpu,
                mem = resource_lumpy_anno.mem
        }
        call variant.concat_vcf as concat_vcf {
            input:
                tools_dir = tools_dir,
                batch_path = batch_path,
                batch_number = batch_number,
                batch_part_vcfs = split_vcfs_flatten,
                batch_part_gvcfs = split_gvcfs_flatten,
                batch_part_vcf_tbis = split_vcf_tbis_flatten,
                batch_part_gvcf_tbis = split_gvcf_tbis_flatten,
                cpu = resource_concat_vcf.cpu,
                mem = resource_concat_vcf.mem
        }
        call variant.variant_mt as variant_mt {
            input:
                tools_dir = tools_dir,
                batch_path = batch_path,
                batch_number = batch_number,
                concat_vcf = concat_vcf.concat_vcf,
                concat_gvcf = concat_vcf.concat_gvcf,
                cpu = resource_variant_mt.cpu,
                mem = resource_variant_mt.mem
        }
        call variant.vqsr_snp as vqsr_snp {
            input:
                tools_dir = tools_dir,
                batch_path = batch_path,
                batch_number = batch_number,
                ref_fa = ref_fa,
                concat_vcf = concat_vcf.concat_vcf,
                concat_vcf_tbi = concat_vcf.concat_vcf_tbi,
                cpu = resource_vqsr_snp.cpu,
                #cpu = 2,
                mem = resource_vqsr_snp.mem
                #mem = 1000
        }
        call variant.vqsr_indel as vqsr_indel {
            input:
                tools_dir = tools_dir,
                batch_path = batch_path,
                batch_number = batch_number,
                ref_fa = ref_fa,
                concat_vcf = concat_vcf.concat_vcf,
                concat_vcf_tbi = concat_vcf.concat_vcf_tbi,
                cpu = resource_vqsr_indel.cpu,
                mem = resource_vqsr_indel.mem
        }
        call variant.concat_snp_indel as concat_snp_indel {
            input:
                tools_dir = tools_dir,
                batch_path = batch_path,
                batch_number = batch_number,
                snp_sort_vcf = vqsr_snp.snp_sort_vcf,
                indel_sort_vcf = vqsr_indel.indel_sort_vcf,
                cpu = resource_concat_snp_indel.cpu,
                mem = resource_concat_snp_indel.mem
        }

        call variant.concat_vcf_2 as concat_vcf_2 {
            input:
                tools_dir = tools_dir,
                batch_path = batch_path,
                batch_number = batch_number,
                batch_part_vcfs_list = concat_vcf.batch_part_vcfs_list,
                batch_part_gvcfs_list = concat_vcf.batch_part_gvcfs_list,
                cpu = resource_concat_vcf_2.cpu,
                mem = resource_concat_vcf_2.mem
        }
        call advance.triploid as triploid {
            input:
                tools_dir = tools_dir,
                conda_path = conda_triploid_env, #此路径可配置在conf文件中
                batch_path = batch_path,
                batch_number = batch_number,
                concat_vcf = concat_vcf_2.concat_vcf,
                concat_vcf_tbi = concat_vcf_2.concat_vcf_tbi,
                cpu = resource_triploid.cpu,
                mem = resource_triploid.mem
        }
        if (! cnvloh_isrun){
            call advance.loh as loh {
                input:
                    concat_vcf = concat_vcf_2.concat_vcf,
                    tools_dir = tools_dir,
                    batch_path = batch_path,
                    batch_number = batch_number,
                    cpu = resource_loh.cpu,
                    mem = resource_loh.mem
            }
        }
        if (cnvloh_isrun){
            call advance.loh_cnv as loh_cnv {
                input:
                    larg_cnv = cnvnator_mops_combine.larg_cnv,
                    concat_vcf = concat_vcf_2.concat_vcf,
                    tools_dir = tools_dir,
                    batch_path = batch_path,
                    batch_number = batch_number,
                    cpu = resource_loh.cpu,
                    mem = resource_loh.mem
            }
        }
        Array[File] bqsr_check_files = split_chr_jobs.bqsr_check_files 
        Array[File] bam_files = split_chr_jobs.bam_files 
        Array[File] bai_files = split_chr_jobs.bai_files 
        Array[Array[String]] split_bed_outs = split_chr_jobs.split_bed_outs 
        call bqsr.all_chr_merge_dup as all_chr_merge_dup {
            input:
                tools_dir = tools_dir,
                bqsr_check_files = bqsr_check_files,
                bam_list = bam_files,
                batch_number = batch_number,
                batch_path = batch_path,
                cpu = resource_all_check_file.cpu,
                mem = resource_all_check_file.mem
        }
        call bqsr.merge_bam_check as merge_bam_check {
            input:
                batch_path = batch_path,
                batch_number = batch_number,
                merge_bam = all_chr_merge_dup.merge_bam,
                cpu = resource_merge_bam_check.cpu,
                mem = resource_merge_bam_check.mem
        }

        call bqsr.qc_collect as qc_collect {
            input:
                batch_path = batch_path,
                batch_number = batch_number,
                tools_dir = tools_dir,
                sex = sex,
                split_bed_outs = split_bed_outs,
                cpu = resource_qc_collect.cpu,
                mem = resource_qc_collect.mem
        }

        Array[Array[File]] sort_bams = qc_align.sort_bams 
        Array[Array[File]] sort_bais = qc_align.sort_bais
        Array[File] sort_bams_flatten = flatten(sort_bams)
        Array[File] sort_bais_flatten = flatten(sort_bais)
        call advance.str_hunter as str_hunter {
            input:
                tools_dir = tools_dir,
                batch_path = batch_path,
                batch_number = batch_number,
                sort_bams = sort_bams_flatten,
                sort_bais = sort_bais_flatten,
                ref_fa = ref_fa,
                sex_txt = qc_collect.sex_txt,
                bam_merge_check = merge_bam_check.bam_merge_check,
                cpu = resource_str_hunter.cpu,
                mem = resource_str_hunter.mem
        }
        call advance.intrauterine_infection as intrauterine_infection {
            input:
                tools_dir = tools_dir,
                bam_merge_check = merge_bam_check.bam_merge_check,
                batch_path = batch_path,
                batch_number = batch_number,
                cpu = resource_intrauterine_infection.cpu,
                mem = resource_intrauterine_infection.mem
        }

        call variant.exon_depth as exon_depth {
            input:
                tools_dir = tools_dir,
                batch_path = batch_path,
                qc_result_txt = qc_collect.qc_result_txt,
                merge_bam = all_chr_merge_dup.merge_bam,
                bam_merge_check = merge_bam_check.bam_merge_check,
                batch_number = batch_number,
                cpu = resource_exon_depth.cpu,
                mem = resource_exon_depth.mem
        }
        call variant.sma as sma {
            input:
                tools_dir = tools_dir,
                bam_list = bam_files,
                bai_list = bai_files,
                batch_path = batch_path,
                cpu = resource_sma.cpu,
                mem = resource_sma.mem
        } 
        scatter (chr in chrs) {
            call variant.cnvnator as cnvnator {
                input:
                    tools_dir = tools_dir,
                    bqsr_check_files = bqsr_check_files,
                    cnvnator_path = cnvnator_path,
                    chr = chr,
                    batch_path = batch_path,
                    cpu = resource_cnvnator.cpu,
                    mem = resource_cnvnator.mem
            }
        }
        scatter (chr in noM_chrs) {
            call variant.mops as mops {
                input:
                    tools_dir = tools_dir,
                    batch_path = batch_path,
                    batch_number = batch_number,
                    chr = chr,
                    bqsr_check_files = bqsr_check_files,
                    qc_result_txt = qc_collect.qc_result_txt,
                    cpu = resource_mops.cpu,
                    mem = resource_mops.mem
            }
        }
        Array[File] all_chr_cnvs = cnvnator.chr_cnv
        Array[File] all_cnmops_cnvs = mops.cnmops_cnv
        call variant.cnvnator_mops_combine as cnvnator_mops_combine {
            input:
                tools_dir = tools_dir,
                cnvnator_results = all_chr_cnvs,
                mops_results = all_cnmops_cnvs,
                batch_path = batch_path,
                batch_number = batch_number,
                cpu = resource_cnvnator_mops_combine.cpu,
                mem = resource_cnvnator_mops_combine.mem
        }
        call report.report as report {
            input:
                str_xlsx = str_hunter.str_xlsx,
                un_rand_bam = un_rand_merge.un_rand_bam,
                unmap_bam = unmap_merge.unmap_bam,
                species_abundance = intrauterine_infection.species_abundance,
                cnv_calls_anno = exon_depth.cnv_calls_anno,
                sma_txt = sma.sma_txt,
                inhouse_wgs_files = inhouse_wgs_files,
                triploid_pdf = triploid.triploid_pdf,
                cnv_loh_tsv = loh_cnv.cnv_loh_tsv,
                loh_tsv = loh.loh_tsv,
                anno_xls = batchCNV_anno.anno_xls,
                anno_result_xlsx = lumpy_anno.anno_result_xlsx,
                anno_xlsx = variant_mt.anno_xlsx,
                database_path = database_path,
                go_path = go_path,
                tools_dir = tools_dir,
                batch_path = batch_path,
                batch_number = batch_number,
                cpu = resource_report.cpu,
                mem = resource_report.mem
        }
    }
}
