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

lines = schema.split(b"\n")
index = 0
for idx, line in enumerate(lines):
    if line.strip().startswith(b'<import'):
        index = idx
        break

lines.insert(index+1, b'<import namespace="http://www.opengis.net/citygml/xAL/2.0" schemaLocation="http://schemas.opengis.net/citygml/xAL/xAL.xsd" />')
schema = b"\n".join(lines)

xsd = etree.fromstring(schema)

xmlschema = etree.XMLSchema(xsd)

doc = etree.parse(os.path.join(SIMULATIONS_DIR,case,case+".gml"))
valid = doc.xmlschema(xsd)

if valid:
    print("The document is valid!")
else:
    print("Invalid document.")
    xmlschema.assert_(doc)