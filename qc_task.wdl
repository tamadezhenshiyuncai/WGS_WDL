# check fastq 这一步会在生成input.json那个外部脚本中进行检测 make_json.py
version development

task filter_split {
    input {
        File read_one
        File read_two
        String tools_dir
        Int cpu
        Float? mem
        String batch_path
        String sample_name
        #Array[String] pieces = ["0001","0002","0003","0004","0005","0006","0007","0008","0009","0010",
        #    "0011","0012","0013","0014","0015","0016","0017","0018","0019","0020"]
    }
    command <<<
        export PATH="~{tools_dir}/tools:$PATH"
        mkdir -p ~{batch_path}/filter/~{sample_name}
        fastp \
          --thread 10 \
          -i ~{read_one} \
          -I ~{read_two} \
          -o ~{batch_path}/filter/~{sample_name}/~{'fastp.' + sample_name}_1.fq.gz \
          -O ~{batch_path}/filter/~{sample_name}/~{'fastp.' + sample_name}_2.fq.gz \
          -s 20 \
          --detect_adapter_for_pe \
          --trim_front1 0 \
          --trim_tail1 0 \
          --trim_front2 0 \
          --trim_tail2 0 \
          -j ~{batch_path}/filter/~{sample_name}/~{sample_name}.fastp.json \
          -h ~{batch_path}/filter/~{sample_name}/~{sample_name}.fastp.html
    >>>
    runtime {
        cpu: cpu
        memory: "~{mem}GB" 
    }
    output {
        Array[Pair[String,Pair[File,File]]] piece_paired_reads = [
            ("0001",("~{batch_path}/filter/~{sample_name}/~{'0001.fastp.' + sample_name + '_1.fq.gz'}","~{batch_path}/filter/~{sample_name}/~{'0001.fastp.' + sample_name + '_2.fq.gz'}")),
            ("0002",("~{batch_path}/filter/~{sample_name}/~{'0002.fastp.' + sample_name + '_1.fq.gz'}","~{batch_path}/filter/~{sample_name}/~{'0002.fastp.' + sample_name + '_2.fq.gz'}")),
            ("0003",("~{batch_path}/filter/~{sample_name}/~{'0003.fastp.' + sample_name + '_1.fq.gz'}","~{batch_path}/filter/~{sample_name}/~{'0003.fastp.' + sample_name + '_2.fq.gz'}")),
            ("0004",("~{batch_path}/filter/~{sample_name}/~{'0004.fastp.' + sample_name + '_1.fq.gz'}","~{batch_path}/filter/~{sample_name}/~{'0004.fastp.' + sample_name + '_2.fq.gz'}")),
            ("0005",("~{batch_path}/filter/~{sample_name}/~{'0005.fastp.' + sample_name + '_1.fq.gz'}","~{batch_path}/filter/~{sample_name}/~{'0005.fastp.' + sample_name + '_2.fq.gz'}")),
            ("0006",("~{batch_path}/filter/~{sample_name}/~{'0006.fastp.' + sample_name + '_1.fq.gz'}","~{batch_path}/filter/~{sample_name}/~{'0006.fastp.' + sample_name + '_2.fq.gz'}")),
            ("0007",("~{batch_path}/filter/~{sample_name}/~{'0007.fastp.' + sample_name + '_1.fq.gz'}","~{batch_path}/filter/~{sample_name}/~{'0007.fastp.' + sample_name + '_2.fq.gz'}")),
            ("0008",("~{batch_path}/filter/~{sample_name}/~{'0008.fastp.' + sample_name + '_1.fq.gz'}","~{batch_path}/filter/~{sample_name}/~{'0008.fastp.' + sample_name + '_2.fq.gz'}")),
            ("0009",("~{batch_path}/filter/~{sample_name}/~{'0009.fastp.' + sample_name + '_1.fq.gz'}","~{batch_path}/filter/~{sample_name}/~{'0009.fastp.' + sample_name + '_2.fq.gz'}")),
            ("0010",("~{batch_path}/filter/~{sample_name}/~{'0010.fastp.' + sample_name + '_1.fq.gz'}","~{batch_path}/filter/~{sample_name}/~{'0010.fastp.' + sample_name + '_2.fq.gz'}")),
            ("0011",("~{batch_path}/filter/~{sample_name}/~{'0011.fastp.' + sample_name + '_1.fq.gz'}","~{batch_path}/filter/~{sample_name}/~{'0011.fastp.' + sample_name + '_2.fq.gz'}")),
            ("0012",("~{batch_path}/filter/~{sample_name}/~{'0012.fastp.' + sample_name + '_1.fq.gz'}","~{batch_path}/filter/~{sample_name}/~{'0012.fastp.' + sample_name + '_2.fq.gz'}")),
            ("0013",("~{batch_path}/filter/~{sample_name}/~{'0013.fastp.' + sample_name + '_1.fq.gz'}","~{batch_path}/filter/~{sample_name}/~{'0013.fastp.' + sample_name + '_2.fq.gz'}")),
            ("0014",("~{batch_path}/filter/~{sample_name}/~{'0014.fastp.' + sample_name + '_1.fq.gz'}","~{batch_path}/filter/~{sample_name}/~{'0014.fastp.' + sample_name + '_2.fq.gz'}")),
            ("0015",("~{batch_path}/filter/~{sample_name}/~{'0015.fastp.' + sample_name + '_1.fq.gz'}","~{batch_path}/filter/~{sample_name}/~{'0015.fastp.' + sample_name + '_2.fq.gz'}")),
            ("0016",("~{batch_path}/filter/~{sample_name}/~{'0016.fastp.' + sample_name + '_1.fq.gz'}","~{batch_path}/filter/~{sample_name}/~{'0016.fastp.' + sample_name + '_2.fq.gz'}")),
            ("0017",("~{batch_path}/filter/~{sample_name}/~{'0017.fastp.' + sample_name + '_1.fq.gz'}","~{batch_path}/filter/~{sample_name}/~{'0017.fastp.' + sample_name + '_2.fq.gz'}")),
            ("0018",("~{batch_path}/filter/~{sample_name}/~{'0018.fastp.' + sample_name + '_1.fq.gz'}","~{batch_path}/filter/~{sample_name}/~{'0018.fastp.' + sample_name + '_2.fq.gz'}")),
            ("0019",("~{batch_path}/filter/~{sample_name}/~{'0019.fastp.' + sample_name + '_1.fq.gz'}","~{batch_path}/filter/~{sample_name}/~{'0019.fastp.' + sample_name + '_2.fq.gz'}")),
            ("0020",("~{batch_path}/filter/~{sample_name}/~{'0020.fastp.' + sample_name + '_1.fq.gz'}","~{batch_path}/filter/~{sample_name}/~{'0020.fastp.' + sample_name + '_2.fq.gz'}")),
        ]
    }
}
