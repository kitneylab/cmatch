#!/usr/bin/env bash

# Create directories for logs and results ouptut
mkdir -p logs;
mkdir -p output_results;

# Activate virtual environment
source ./venv/bin/activate;

# Launch test for algorithm CM_0
#python Testing_Algorithm_CM_0.py;

# Launch tests (3 loops) for CM_1 with Vio-0000 of length x2 and x3
#python cm1_cat2.py;
#python cm1_cat2min.py;
#python cm1_cat3.py;
#python cm1_cat3min.py;

# Launch tests (10 loops) for CM_2 with Vio-0000 of length x2 and x3
#python cm2_cat2.py;
#python cm2_cat2min.py;
#python cm2_cat3.py;
#python cm2_cat3min.py;

python cm1_ten.py
