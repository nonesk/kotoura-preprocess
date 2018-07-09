#! /usr/bin/env python3

from libsbml import *
import pandas as pd
import sys
import re

from collections import defaultdict, OrderedDict

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

def getSpecies(model, typeIds):
    species = dict()
    for sp in model.getListOfSpecies():
        if sp.getSpeciesType() in typeIds and sp.getCompartment() == "c_03":
            species[sp.getSpeciesType()] = sp.getId()
    return species


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

ids = OrderedDict()
missing = set(targetNames)
for sp in species :
    for target in targetNames:
        if re.match("^"+re.escape(sp.getName())+"$",
                    translate(target, translations),
                    re.IGNORECASE):
            print(sp.getName(), sp.getId())
            ids[target] = sp.getId()
            if target in missing :
                missing.remove(target)

print(ids)
print(len(ids))

print(missing)

species = getSpecies(model, ids.values())


# print(set([sp.split()[-1] for sp in sum(species.values())]))

print (species)

print(len(species))

sIds = [species[ids[target]] for target in targetNames]
print(sIds)

se = pd.Series(sIds)
df['meta', 'sid'] = se
print(df)
df = df.sort_index(1, 0, ascending=False)
df.columns = ['_'.join(col) for col in df.columns.values]
df.to_csv("processed_" + csv)

print(model.getNumSpecies())
print(model.getNumReactions())
print(model.getNumCompartments())
