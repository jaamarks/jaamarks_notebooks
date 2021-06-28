#!/usr/bin/bash

base_dir=$PWD
method_name=_methods-fwgwas-naive-scz2-20180911.py

~/bin/qsub_job.sh \
    --job_name naive-pipeline \
    --script_prefix  $base_dir/naive-scz2-gw \
    --mem 10 \
    --cpu 8 \
    --priority 0 \
    --program python $base_dir/$method_name
