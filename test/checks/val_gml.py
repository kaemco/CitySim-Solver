#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Validate GML file using XSD schema.

Created on Fri Feb 10 15:00:00 2023

@author: gperonato
"""

from lxml import etree
import os
import requests

SIMULATIONS_DIR = "../simulations"
case = "600_profiles"

schema = requests.get("http://www.sig3d.org/citygml/2.0/energy/1.0/EnergyADE.xsd").content

doc = etree.parse(os.path.join(SIMULATIONS_DIR,case,case+".gml"))

xsd = etree.fromstring(schema)

xmlschema = etree.XMLSchema(xsd)
valid = doc.xmlschema(xsd)

if valid:
    print("The document is valid!")
else:
    print("Invalid document.")
    xmlschema.assert_(doc)

