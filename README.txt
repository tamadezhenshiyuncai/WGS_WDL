最终外面会写一个脚本去生成input.json文件，这个脚本会去解析 family3.fq.ls 这个fastq的配置文件
task check_fq 只要进行单个的数据存在与否判断就可以
在没有容器前，需要一个conf文件记录一些存放配置文件的路径 如： tools_dir 路径 GO路径等
