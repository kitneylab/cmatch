#!/usr/bin/env bash

source ./venv/bin/activate;
#python run_algo0.py;
#python run_algo1.py;
#python run_algo2.py;
#python run_algo1_lycopene.py
#python reconstruction.py

# 1 vs 1
#python test_cm2_vio_easy_1vs1_th99.py
#python test_cm2_vio_hard_1vs1_th75.py
#python test_cm2_vio_hard_1vs1_th99.py
#python test_cm2_vio_easy_1vs1_th75.py

# 1 vs ALL
#python test_cm2_vio_easy_1vsAll_th99.py
#python test_cm2_vio_easy_1vsAll_th75.py
#python test_cm2_vio_hard_1vsAll_th99.py
#python test_cm2_vio_hard_1vsAll_th75.py

# Lycopene

#python test_cm2_lycopene_1.py
python test_cm2_lycopene_sanger_100.py
