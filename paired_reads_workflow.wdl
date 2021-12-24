version development

import "structs.wdl" as structs
    alias Modules as Modules
    alias Resource as Resource
    alias FileFilePair as FileFilePair
    alias StringArrayFilePair as StringArrayFilePair
    alias StringArrayStringPair as StringArrayStringPair
import "qc_task.wdl" as qc
import "align_piece_workflow.wdl" as align_piece_workflow

workflow paired_reads_fastp {
    input {
        String batch_path
        String batch_number
        String tools_dir
        String project_path
        String ref_fa
        Int qc_cpu
        Float? qc_mem
        Int align_cpu
        Float? align_mem
        Int sort_cpu
        Float? sort_mem
        #Array[Pair[File,File]] batch_paired_reads  
        Array[FileFilePair] batch_paired_reads  
    }
    scatter (pair_reads in batch_paired_reads) {
        File read_one = pair_reads.left
        File read_two = pair_reads.right
        String read_one_base = basename(read_one,"_1.fq.gz")
        String sample_name = sub(read_one_base,"fastp.","")
        call qc.filter_split as filter_split {
            input:
                read_one = read_one,
                read_two = read_two,
                tools_dir = tools_dir,
                cpu = qc_cpu,
                mem = qc_mem,
                batch_path = batch_path,
                sample_name = sample_name
        }
        call align_piece_workflow.align_piece as align_piece {
            input:
                tools_dir = tools_dir,
                project_path = project_path,
                batch_number = batch_number,
                batch_path = batch_path,
                ref_fa = ref_fa,
                sample_name = sample_name,
                piece_paired_reads = filter_split.piece_paired_reads,
                align_cpu = align_cpu,
                align_mem = align_mem,
                sort_cpu = sort_cpu,
                sort_mem = sort_mem
        }
        String align_path = "~{batch_path}/align/~{sample_name}"
    }
    output {
        Array[String] align_path_arr = align_path
        Array[Array[String]] done_flag = align_piece.done_flag
        Array[Array[File]] sort_bams = align_piece.sort_bams
        Array[Array[File]] sort_bais = align_piece.sort_bais
    }
}
