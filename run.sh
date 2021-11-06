#! /user/bin/env bash

# Create directories for log and results ouptut
mkdir logs;
mkdir output_results;

# Activate virtual environment
source ./venv/bin/activate;

# Launch test for algorithm CM_0
python Testing_Algorithm_CM_0.py;
