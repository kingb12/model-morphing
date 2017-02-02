kbase-login bking -p sea:krhsUW
python ./lib/model-morphing/scripts/prepare_deploy_cfg.py ./deploy.cfg ./lib/model-morphing/kbconfig.properties
python ./lib/model-morphing/scripts/extract_token.py ~/.kbase_config ./work/token
export KB_AUTH_TOKEN=`cat /kb/module/work/token`

python ./lib/model-morphing/mari_to_janna.py $1 &
python ./lib/model-morphing/mari_to_bark.py $2
