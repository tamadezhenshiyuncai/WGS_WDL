#! /usr/bin/env python
#-*-coding:utf-8-*-
__version__ = 'v2'
__date__ = '2021/12/06'
#Modify by fuxiangke 第二版主要修改了一个struct变量中的结构，因为染色体切份有按6M切份的，有按24M切份的

import os
import json
import argparse
from configparser import ConfigParser, ExtendedInterpolation
from enum import Enum
from collections import defaultdict

class CPU(Enum):
    cnvnator_mops_combine = 10
    variant_mt = 7
    align_sort = 1
    loh = 1
    mops = 2
    bqsr_check = 1
    all_chr_merge_dup = 2
    unmap_merge = 5
    batchCNV = 10
    qc_collect = 10
    align = 6
    concat_vcf_2 = 1
    vqsr_indel = 2
    merge_bam_check = 1
    batchCNV_anno = 1
    exon_depth = 1
    filter_split = 4
    concat_vcf = 1
    lumpy_chr_disco = 1
    variant_annotation = 7
    cnvnator = 2
    qc_bed_split = 1
    bam_chr_part_split = 1
    sma = 1
    report = 2
    bqsr = 1
    final_lumpy = 1
    un_rand_merge = 5
    lumpy_anno = 20
    variant_hc = 2
    concat_snp_indel = 1
    batchCNV_breakpoint = 1
    intrauterine_infection = 5
    str_hunter = 4
    triploid = 6
    lumpy_chr_split = 1
    bam_merge_dup = 2
    vqsr_snp = 2

class MEM(Enum):
    cnvnator_mops_combine = 15.0
    variant_mt = 5.0
    align_sort = 3.0
    loh = 3.0
    mops = 40.0
    bqsr_check = 1.0
    all_chr_merge_dup = 30.0
    unmap_merge = 7.0
    batchCNV = 24.0
    qc_collect = 2.0
    align = 6.0
    concat_vcf_2 = 3.0
    vqsr_indel = 6.0
    merge_bam_check = 1.0
    batchCNV_anno = 3.0
    exon_depth = 3.0
    filter_split = 1.0
    concat_vcf = 4.0
    lumpy_chr_disco = 3.0
    variant_annotation = 13.0
    cnvnator = 4.0
    qc_bed_split = 3.0
    bam_chr_part_split = 3.0
    sma = 10.0
    report = 3.0
    bqsr = 10.0
    final_lumpy = 3.0
    un_rand_merge = 7.0
    lumpy_anno = 11.0
    variant_hc = 6.0
    concat_snp_indel = 6.0
    batchCNV_breakpoint = 3.0
    intrauterine_infection = 5.0
    str_hunter = 8.0
    triploid = 8.0
    lumpy_chr_split = 3.0
    bam_merge_dup = 40.0
    vqsr_snp = 9.0

def init_args():
    parser = argparse.ArgumentParser(description='读取sample.lst,与配置文件config.ini,生成WDL流程的输入文件input.json')
    parser.add_argument('-l', '--sample_list', type=str, required=True, help='输入的sample.lst文件,该文件包含六列信息')
    parser.add_argument('-c', '--config', type=str, required=True, help='配置文件所在的路径,包含了input.json文件中某些输入项的值')
    args = parser.parse_args()
    return args

def _isexists(fastq):
    """检查文件是否存在，不存在则抛出IO错误"""
    if not os.path.exists(fastq):
        raise IOError("{fastq} is not exists".format(fastq=fastq))

def parse_smp_file(sample_list):
    """sample_list文件的fastq路径只填写read1的路径,通过read1在同路径下找read2
    另外,虽然WDL流程本身是会根据类型去check数据存在与否，此处还是再进行check一次
    """
    reads_info_list = []
    reads_info_dict = defaultdict(dict)
    with open(sample_list, "r") as inf:
        for each in inf:
            cols = each.strip().split()
            sample_name, read1_fastq = cols[:2]
            read2_fastq = read1_fastq.replace("_1","_2")
            _isexists(read1_fastq)
            _isexists(read2_fastq)
            _dict = reads_info_dict[sample_name]
            if not _dict.get("left"):
                _dict["left"] = sample_name
            if not _dict.get("right"):
                _dict["right"] = [{"left":read1_fastq,"right":read2_fastq}]
            else:
                _dict["right"].append({"left":read1_fastq,"right":read2_fastq})
    for each in reads_info_dict:
        reads_info_list.append(reads_info_dict[each])
    return reads_info_list

def parse_config_file(conf_file):
    """从配置的ini文件中读取流程所需要的信息,有任何改动也可在这个文件中修改"""
    cfg = ConfigParser(interpolation=ExtendedInterpolation())
    cfg.read(conf_file)
    conf_dict = {}
    conf_dict["wgs_pip.tools_dir"] = cfg.get('Tools', 'tools_dir')
    conf_dict["wgs_pip.sex"] = cfg.get('Project', 'sex')
    conf_dict["wgs_pip.go_path"] = cfg.get('Tools', 'go_path')
    conf_dict["wgs_pip.conda_triploid_env"] = cfg.get('Conda', 'triploid')
    conf_dict["wgs_pip.indels_sites_vcf"] = cfg.get('Vcfs', 'indel_vcf')
    conf_dict["wgs_pip.project_path"] = cfg.get('Project', 'project_path')
    conf_dict["wgs_pip.conda_lumpy_env"] = cfg.get('Conda', 'lumpy')
    conf_dict["wgs_pip.dbsnp_vcf"] = cfg.get('Vcfs', 'dbsnp_vcf')
    conf_dict["wgs_pip.high_confidence_sites_vcf"] = cfg.get('Vcfs', 'confidence_vcf')
    conf_dict["wgs_pip.cnvnator_path"] = cfg.get('Tools', 'cnvnator_path')
    conf_dict["wgs_pip.omni_sites_vcf"] = cfg.get('Vcfs', 'omni_vcf')
    conf_dict["wgs_pip.ref_fa"] = cfg.get('Project', 'ref_fa')
    conf_dict["wgs_pip.database_path"] = cfg.get('Project', 'database_path')
    return cfg, conf_dict

def chr_split_beds(bed6m_dir, bed24m_dir):
    """拆分染色体"""
    chr_parts_list = []
    chr_dict = defaultdict(dict)
    for each in os.listdir(bed6m_dir):
        if each.endswith(".bed"):
            chrom, _, _ = each.rsplit('.',2)
            bed_abspath = os.path.join(bed6m_dir,each)
            if not chr_dict[chrom].get("left"):
                chr_dict[chrom]["left"] = chrom
            if not chr_dict[chrom].get("right"):
                chr_dict[chrom]["right"] = [bed_abspath]
            else:
                chr_dict[chrom]["right"].append(bed_abspath)
    for each in os.listdir(bed24m_dir):
        if each.endswith(".bed"):
            chrom, _, _ = each.rsplit('.',2)
            bed_abspath = os.path.join(bed24m_dir,each)
            if not chr_dict[chrom].get("left"):
                chr_dict[chrom]["left"] = chrom
            if not chr_dict[chrom].get("mid"):
                chr_dict[chrom]["mid"] = [bed_abspath]
            else:
                chr_dict[chrom]["mid"].append(bed_abspath)
    for each_chr in chr_dict:
        chr_parts_list.append(chr_dict[each_chr])
    return chr_parts_list

def main():
    json_dict = {}
    _dict = {}
    json_dict["wgs_pip.resources"]  = _dict
    for each in zip(CPU.__members__.items(),MEM.__members__.items()):
        name = each[0][0] 
        _dict.update({name: {"cpu": each[0][1].value,"mem": each[1][1].value}})
    args = init_args()
    conf_file = args.config
    sample_list = args.sample_list
    cfg, conf_dict = parse_config_file(conf_file)
    chr_split6m_path = cfg.get('Project', 'chr_parts6m_path')
    chr_split24m_path = cfg.get('Project', 'chr_parts24m_path')
    chr_parts_list = chr_split_beds(chr_split6m_path,chr_split24m_path)
    json_dict["wgs_pip.noM_chrs"] = ["chr1","chr2","chr3","chr4","chr5","chr6","chr7",
    "chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18",
    "chr19","chr20","chr21","chr22","chrX","chrY"]
    json_dict["wgs_pip.cnvloh_isrun"] = False
    json_dict.update(conf_dict)
    json_dict["wgs_pip.chr_parts"] = chr_parts_list
    reads_info_list = parse_smp_file(sample_list)
    json_dict["wgs_pip.reads_info"] = reads_info_list
    json_dict["wgs_pip.chrs"] = ["chr1","chr2","chr3","chr4","chr5","chr6","chr7",
    "chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17",
    "chr18","chr19","chr20","chr21","chr22","chrX","chrY","chrM_NC_012920.1"]
    return json_dict

if __name__ == "__main__":
    json_dict = main()
    with open("input.json", "w") as outf:
        outf.write(json.dumps(json_dict, indent=4))
