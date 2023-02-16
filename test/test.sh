#!/bin/bash

echo 'Running simulations...'
cd simulations
bash ./run.sh
echo ''
cd ../checks
python3 check.py
echo 'Validating GML...'
python3 val_gml.py
