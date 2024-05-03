#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 27 18:28:00 2022

@author: gperonato
"""

import pandas as pd
import os

SIMULATIONS_DIR = "../simulations"
case = "600"

# This is the case we want to test
results = pd.read_csv(os.path.join(SIMULATIONS_DIR,case,"{}_TH.out".format(case)),
                      sep='\s+',index_col=0)

# These are pre-computed reference results
reference = pd.read_csv(os.path.join("reference","{}_TH.out".format(case)),
                      sep='\s+',index_col=0)


# These are reference results pre-computed from BESTEST
besttest = pd.read_csv(os.path.join("reference","BESTEST.csv"),index_col=0)


# Load base results
heating = results.loc[:,results.columns.str.contains("Heating")]/1000
cooling = results.loc[:,results.columns.str.contains("Cooling")].abs()/1000


# Checks    
print("Checking results...")
if not besttest[case]["annual_heating_min"] < heating.sum().values[0]/1000 < besttest[case]["annual_heating_max"]:
    raise ValueError("Annual heating is outside the BESTTEST range")
    
if not besttest[case]["annual_cooling_min"] < cooling.sum().values[0]/1000 < besttest[case]["annual_cooling_max"]:
    raise ValueError("Annual cooling is outside the BESTTEST range")
   
if not besttest[case]["peak_heating_min"] < heating.max().values[0] < besttest[case]["peak_heating_max"]:
    raise ValueError("Peak heating is outside the BESTTEST range")
    
if not besttest[case]["peak_cooling_min"] < cooling.max().values[0] < besttest[case]["peak_cooling_max"]:
    raise ValueError("Peak cooling is outside the BESTTEST range")
    
pd.testing.assert_frame_equal(reference, results, atol=0.002)
    
print("All checks passed!")
    
