#/zfsyt1/B2C_RD_P2/USER/fuxiangke/software/Miniconda2/envs/cromwell/bin/cromwell submit \
#        -i /zfsyt1/B2C_RD_P2/USER/fuxiangke/wgs_server_mode_1222_v2/input.json \
#        -o /zfsyt1/B2C_RD_P2/USER/fuxiangke/wgs_server_mode_1222_v2/options.json \
#        -h http://192.168.3.6:5002 \
#        -p /zfsyt1/B2C_RD_P2/USER/fuxiangke/WGS_WDL_SGE_V3/wgs_tasks.zip \
#        /zfsyt1/B2C_RD_P2/USER/fuxiangke/WGS_WDL_SGE_V3/wgs_workflow_v3.wdl 
/zfsyt1/B2C_RD_P2/USER/fuxiangke/software/Miniconda2/envs/cromwell/bin/java -jar /zfsyt1/B2C_RD_P2/USER/fuxiangke/software/Miniconda2/envs/cromwell/share/cromwell/cromwell-72.jar  submit \
        -i /zfsyt1/B2C_RD_P2/USER/fuxiangke/wgs_server_mode_1222/input.json \
        -o /zfsyt1/B2C_RD_P2/USER/fuxiangke/wgs_server_mode_1222/options.json \
        -h http://192.168.3.6:5002 \
        -p /zfsyt1/B2C_RD_P2/USER/fuxiangke/WGS_WDL_SGE_V3/wgs_tasks.zip \
        /zfsyt1/B2C_RD_P2/USER/fuxiangke/WGS_WDL_SGE_V3/wgs_workflow_v3.wdl 
