#!/bin/bash

cd simulations
bash ./run.sh
cd ../checks
python3 check.py

