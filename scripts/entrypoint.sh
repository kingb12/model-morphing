#!/bin/bash

# SDKBase set up, not sure what all is needed
. /kb/deployment/user-env.sh

export PYTHONPATH="${PYTHONPATH}:/kb/module/lib/model-morphing/lib/"
export PERL5LIB="/kb/module/lib:${PERL5LIB}:${PATH}:/kb/module/lib/model-morphing/:kb/runtime/bin:/kb/runtime/java/bin"
export KB_DEPLOYMENT_CONFIG="/kb/module/deploy.cfg"

# Login to KBase in this enviroment (creates a .kbase_config file in the home directory)
if [ $# -gt 1 ] && [ "${1}" = "login" ] ; then
  sh ./lib/model-morphing/scripts/login.sh "${@:2}"
else
  echo "Please enter your KBase username (sign up at http://kbase.us):"
  read username
  sh ./lib/model-morphing/scripts/login.sh $username
fi

# Set Up workdir

python ./lib/model-morphing/scripts/prepare_deploy_cfg.py ./deploy.cfg ./lib/model-morphing/kbconfig.properties
python ./lib/model-morphing/scripts/extract_token.py ~/.kbase_config ./work/token
export KB_AUTH_TOKEN=`cat /kb/module/work/token`

if [ $# -eq 0 ] ; then
  bash
elif [ "${1}" = "server" ] ; then
  echo "Run Server"
  sh ./scripts/start-server.sh
  ./scripts/st
elif [ "${1}" = "test" ] ; then
  echo "Run Tests"
  make testr
elif [ "${1}" = "async" ] ; then
  sh ./scripts/run_async.sh
elif [ "${1}" = "init" ] ; then
  echo "Initialize module"
elif [ "${1}" = "bash" ] ; then
  bash
elif [ "${1}" = "report" ] ; then
  export "KB_SDK_COMPILE_REPORT_FILE=./work/compile_report.json"
  make compile
elif [ "${1}" = "login" ] ; then
  ipython -i /kb/module/lib/model-morphing/scripts/._mm_init.py
else
  echo Unknown
fi