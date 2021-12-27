# WGS_WDL

## 流程准备与运行
    1. 用WDL语法编写的workflow
    2. 安装cromwell, cromwell是启动WDL流程的引擎，有两种模式，一种是run模式，一种是server模式，测试小的流程一般用run，真正生产使用为了支持断点续跑还是要选择server模式
    3. 安装mariaDB 本地数据库，安装并启动服务在后台，并将其配置到cromwell的配置conf文件中的database块中，注意：cromwell并不会为我们创建数据库，我们在安装mariaDB之后要自行启动数据库实例，并登入账号，第一个登入的账号即为管理账号; 之后执行如下sql语句：
        CREATE DATABASE IF NOT EXISTS cromwell default charset utf8mb4 COLLATE utf8mb4_unicode_ci;
        CREATE USER 'cromwell'@'localhost' IDENTIFIED BY 'cromwell';
        GRANT all privileges ON cromwell.* TO 'cromwell'@'localhost';
    4. womtool 检查wdl流程的语法
    5. 语法没有问题后，用womtool 生成input.json模板，并用实际数据配置，真正应用到生产时，input.json将会有专门的脚本生成
    6. 将cromwell 服务脚本挂在后台：
nohup /share/nastj10/B2C_RD_P2/USER/fuxiangke/software/Miniconda3/envs/cromwell/bin/cromwell server -Dconfig.file=/share/nastj10/B2C_RD_P2/USER/fuxiangke/software/pipline_wgs/WGS_WDL_SGE/SGE_SERVER.conf &
    7. 执行如下脚本，将任务往服务节点上推送：
`/share/nastj10/B2C_RD_P2/USER/fuxiangke/software/Miniconda3/envs/cromwell/bin/cromwell submit \
    -i /share/nastj4/B2C_RD_P2/USR/fuxiangke/wgs_server_test/new_input.json \
    -o /share/nastj4/B2C_RD_P2/USR/fuxiangke/wgs_server_test/options.json \
    -h http://192.168.136.114:5002 \
    -p /share/nastj10/B2C_RD_P2/USER/fuxiangke/software/pipline_wgs/WGS_WDL_SGE/wgs_tasks.zip \
    /share/nastj10/B2C_RD_P2/USER/fuxiangke/software/pipline_wgs/WGS_WDL_SGE/wgs_workflow.wdl`
<font color=#00ffff size=3>null</font>
