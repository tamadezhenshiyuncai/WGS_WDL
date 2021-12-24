version development

struct Resource {
    Int cpu
    Float? mem
}

struct FileFilePair {
    File left
    File right
}

struct StringArrayFilePair {
    String left
    Array[FileFilePair] right
}

struct StringArrayStringPair {
    String left
    Array[String] mid
    Array[String] right
}

struct Modules {
    Resource filter_split
    Resource align
    Resource align_sort
    Resource un_rand_merge
    Resource unmap_merge
    Resource bam_merge_dup
    Resource bqsr
    Resource bqsr_check
    Resource all_chr_merge_dup
    Resource merge_bam_check
    Resource bam_chr_part_split
    Resource qc_bed_split
    Resource qc_collect
    Resource variant_hc
    Resource concat_vcf
    Resource variant_mt
    Resource vqsr_snp
    Resource vqsr_indel
    Resource concat_snp_indel
    Resource concat_vcf_2
    Resource variant_annotation
    Resource exon_depth
    Resource sma
    Resource cnvnator
    Resource mops
    Resource cnvnator_mops_combine
    Resource lumpy_chr_split
    Resource lumpy_chr_disco
    Resource final_lumpy
    Resource lumpy_anno
    Resource batchCNV
    Resource batchCNV_breakpoint
    Resource batchCNV_anno
    Resource triploid
    Resource loh
    Resource str_hunter
    Resource intrauterine_infection
    Resource report
}
