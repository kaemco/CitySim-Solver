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
missing_import = etree.XML(b'<import xmlns="http://www.w3.org/2001/XMLSchema" xmlns:energy="http://www.sig3d.org/citygml/2.0/energy/1.0" xmlns:app="http://www.opengis.net/citygml/appearance/2.0" xmlns:grp="http://www.opengis.net/citygml/cityobjectgroup/2.0" xmlns:core="http://www.opengis.net/citygml/2.0" xmlns:bldg="http://www.opengis.net/citygml/building/2.0" xmlns:gml="http://www.opengis.net/gml" namespace="http://www.opengis.net/citygml/xAL/2.0" schemaLocation="http://schemas.opengis.net/citygml/xAL/xAL.xsd"/>')

schema = requests.get("http://www.sig3d.org/citygml/2.0/energy/1.0/EnergyADE.xsd").content

xsd = etree.fromstring(schema)
xsd.insert(0, missing_import)
xmlschema = etree.XMLSchema(xsd)

doc = etree.parse(os.path.join(SIMULATIONS_DIR,case,case+".gml"))

if xmlschema.validate(doc):
    print("The document is valid!")
else:
    print("Invalid document.")
    xmlschema.assert_(doc)
