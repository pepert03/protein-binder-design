#!/bin/bash
cd ~/glp1r
. ~/tnf/venv/bin/activate
SP=~/tnf/venv/lib/python3.10/site-packages/nvidia
export LD_LIBRARY_PATH=$SP/cusparse/lib:$SP/cublas/lib:$SP/cudnn/lib:$SP/cufft/lib:$SP/cusolver/lib:$SP/cuda_runtime/lib:$SP/nvjitlink/lib:$LD_LIBRARY_PATH
exec python -u scripts/run_glp1r_colabdesign.py
