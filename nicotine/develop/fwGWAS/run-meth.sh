#!/usr/bin/bash

base_dir=$PWD
method_name=_methods_fwGWAS_20190111.py
~/bin/qsub_job.sh \
    --job_name naive-pipeline \
    --script_prefix  $base_dir/naive-gw \
    --mem 10 \
    --cpu 8 \
    --priority 0 \
    --program python $base_dir/$method_name
