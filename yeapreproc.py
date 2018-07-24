#! /usr/bin/env python3

from libsbml import *
import pandas as pd
import sys
import re

from collections import defaultdict, OrderedDict

from math import sqrt

# Replaces compounds names from CSV to match network names
translations = {
    "Red_Glutathione": "Glutathione",
    "Oxid_Glutathione": "Glutathione disulfide",
    "UDP-Glucose": "UDP-D-Glucose",
    "Isopropylmalate": "2-isopropylmalate",
    "L-aspartic acid": "L-aspartate",
    "Citric acid": "citrate(3-)",
    "Fumaric acid": "fumarate",
    "Malic acid": "(S)-malate",
    "L-glutamic acid": "L-glutamate",
    "Succinic acid": "succinate",
    "L-ornithine": "ornithine",
    "Orotidine": "orotidine 5'-(dihydrogen phosphate)"
}


def compute_delta(df, control, condition):
    """
    Compute concentration variations between control and ethanol
        :param df: pandas dataframe
        :param control: control column name
        :param condition: condition column name
    """
    ctrl = df[control]
    eth = df[condition]
    # avg difference standard deviation at 95%
    deltasd = (eth['sd'].pow(2) / 6 + ctrl['sd'].pow(2) / 6).pow(1 / 2)
    # compute delta min, delta max
    df['delta', 'min'] = eth['avg'] - ctrl['avg'] - 1.96 * deltasd
    df['delta', 'max'] = eth['avg'] - ctrl['avg'] + 1.96 * deltasd
    return df


def fix_headers(df):
    """
    Update dataframe header to allow multiple row header
        :param df: pandas dataframe
    """
    a = df.columns.get_level_values(0).to_series()
    b = a.mask(a.str.startswith("Unnamed")).ffill()
    b[0] = "meta"
    df.columns = pd.MultiIndex.from_arrays([b, df.columns.get_level_values(1)])
    return df


def translate(name, d):
    """
    Replace compound name from dictionnary, if necessary
        :param name: compound name
        :param d: translation dictionnary
    """
    return d.get(name) or name


def getSpecies(model, typeIds):
    """
    Creates dictionnary of species as speciesType => speciesId
        :param model: SBML model
        :param typeIds: list of species type IDs
    """
    species = dict()
    for sp in model.getListOfSpecies():
        if sp.getSpeciesType() in typeIds and sp.getCompartment() == "c_03":
            species[sp.getSpeciesType()] = sp.getId()
    return species


assert len(sys.argv) > 2, "USAGE : yeastpreproc.py <CSV file> <condition column name>"

csv = sys.argv[1]
condition = sys.argv[2]

# Read CSV data
df = pd.read_csv(csv, sep=',',
                 decimal=".",
                 float_precision="high",
                 header=[0, 1])
fix_headers(df)
compute_delta(df, '0perc', condition)

# Read SBML file
reader = SBMLReader()
document = reader.readSBML("yeast_5.01_model.xml")
if document.getNumErrors():
    raise ValueError("SBML File corrupted\n")
model = document.getModel()
species = model.getListOfSpeciesTypes()
# Names of metabolites
targetNames = df["meta", "Metabolites"].tolist()

# get metabolite type IDs and store them as name => SBML type ID
ids = OrderedDict()
missing = set(targetNames)
for sp in species:
    for target in targetNames:
        if re.match("^" + re.escape(sp.getName()) + "$",
                    translate(target, translations),
                    re.IGNORECASE):
            #print(sp.getName(), sp.getId())
            ids[target] = sp.getId()
            if target in missing:
                missing.remove(target)

if missing: 
    raise ValueError("Missing compounds : \n{}".format(missing))

# get species ids
species = getSpecies(model, ids.values())
# get list of species IDs ordered like CSV rows
sIds= [species[ids[target]] for target in targetNames]
# add species IDs to dataframe
se = pd.Series(sIds)
df['meta', 'sid'] = se
# organize dataframe column headers
df = df.sort_index(1, 0, ascending=False)
df.columns = ['_'.join(col) for col in df.columns.values]
print(df)
df.to_csv("processed_" + csv)

