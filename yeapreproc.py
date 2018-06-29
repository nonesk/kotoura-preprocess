#! /usr/bin/env python3

from libsbml import *
import pandas as pd
import sys
import re

from math import sqrt

translations = {
    "Red_Glutathione" : "Glutathione",
    "Oxid_Glutathione" : "Glutathione disulfide",
    "UDP-Glucose" : "UDP-D-Glucose",
    "Isopropylmalate" : "2-isopropylmalate",
    "L-aspartic acid" : "L-aspartate",
    "Citric acid" : "citrate(3-)",
    "Fumaric acid" : "fumarate",
    "Malic acid" : "(S)-malate",
    "L-glutamic acid" : "L-glutamate",
    "Succinic acid" : "succinate",
    "L-ornithine" : "ornithine",
    "Orotidine" : "orotidine 5'-(dihydrogen phosphate)"
}

def compute_delta(df, condition1, condition2):
    ctrl = df[condition1]
    eth = df[condition2]
    # avg difference standard deviation at 95%
    deltasd = (eth['sd'].pow(2)/6 + ctrl['sd'].pow(2)/6).pow(1/2)
    # compute delta min, delta max
    df['delta', 'min'] = eth['avg'] - ctrl['avg'] - 1.96 * deltasd
    df['delta', 'max'] = eth['avg'] - ctrl['avg'] + 1.96 * deltasd
    return df

def fix_headers(df):
    a = df.columns.get_level_values(0).to_series()
    b = a.mask(a.str.startswith("Unnamed")).ffill()
    b[0] = "meta"
    df.columns = pd.MultiIndex.from_arrays([b,df.columns.get_level_values(1)])
    return df

def translate(name, d):
    return d.get(name) or name


csv = sys.argv[1]

df = pd.read_csv(csv, sep=',', decimal=".", float_precision="high", header=[0,1])
fix_headers(df)
compute_delta(df, '0perc', '2perc')
print(df)


reader = SBMLReader()
document = reader.readSBML("yeast_5.01_model.xml")
if document.getNumErrors() :
    raise ValueError("SBML File corrupted\n")
model = document.getModel()
species = model.getListOfSpeciesTypes()
targetNames = df["meta", "Metabolites"].tolist()

ids = set()
missing = set(targetNames)
for sp in species :
    for target in targetNames:
        if re.match("^"+re.escape(sp.getName())+"$",
                    translate(target, translations),
                    re.IGNORECASE):
            print(sp.getName(), sp.getId())
            ids.add(sp.getId())
            if target in missing :
                missing.remove(target)

print(ids)
print(len(ids))

print(missing)

model.getListOfSpeciesTypes()
