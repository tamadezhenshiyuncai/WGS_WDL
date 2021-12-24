version development

import "align_task.wdl" as align

workflow align_piece {
    input {
        String tools_dir
        String project_path
        String batch_number
        String batch_path
        String ref_fa
        String sample_name
        Array[Pair[String,Pair[File,File]]] piece_paired_reads
        Int align_cpu
        Float? align_mem
        Int sort_cpu
        Float? sort_mem
    }
    scatter (each_paired_reads in piece_paired_reads) {
        String piece = each_paired_reads.left
        Pair[File,File] paired_reads = each_paired_reads.right
        call align.align as align_t {
            input:
                tools_dir = tools_dir,
                clean_read_one = paired_reads.left,
                clean_read_two = paired_reads.right,
                batch_number = batch_number,
                split_ref_fa = ref_fa,
                project_path = project_path,
                sample_name = sample_name,
                piece = piece,
                cpu = align_cpu,
                mem = align_mem
        }
        call align.align_sort as align_sort {
            input:
                tools_dir = tools_dir,
                batch_number = batch_number,
                sample_name = sample_name,
                batch_path = batch_path,
                piece = piece,
                piece_out_bam = align_t.out_bam,
                cpu = sort_cpu,
                mem = sort_mem
        }
    }
    output {
        Array[String] done_flag = align_sort.done_flag
        Array[File] sort_bams = align_sort.sort_bam
        Array[File] sort_bais = align_sort.sort_bai
    }
}
