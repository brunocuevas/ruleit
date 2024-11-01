import streamlit as st
import pandas as pd
import json
from rdkit.Chem import rdChemReactions as rdr
from ruleit.hash import hash_molecules
import numpy as np


def generate_stochiometry_dict(reaction):
    r = rdr.ReactionFromSmarts(reaction, useSmiles=True)
    reactants = r.GetReactants()
    products = r.GetProducts()
    s = dict()
    for mol in reactants:
        try:
            s[hash_molecules(mol=mol)] -= 1
        except KeyError:
            s[hash_molecules(mol=mol)] = -1

    for mol in products:
        try:
            s[hash_molecules(mol=mol)] += 1
        except KeyError:
            s[hash_molecules(mol=mol)] = 1
    return s


st.title("Rule-It: Stoichiometry 😳")

"""
This page allows to convert the outputs of the other pages (´prune´ and ´expand´) into the input of
an autocatalytic-cycle detection tool, freely available at [https://github.com/vblancoOR/autocatatalyticsubnetworks](https://github.com/vblancoOR/autocatatalyticsubnetworks)
"""

upload = st.file_uploader(label="upload")
upload_name = upload.name
run = st.button('Run')
flag_too_big = False
if upload and run:
    name = upload.name
    upload = json.load(upload)
    reactions = upload['reactions']
    molecules = upload['seeds']
    molecules_index = list(map(hash_molecules, molecules))
    rules = upload['rules']
    st.metric('#reactions', len(reactions))

    if len(reactions) > 1000:
        st.write("the stoichiometric matrix is too big to display")
        flag_too_big = True

    stoichiometries = []
    for r in reactions:
        try:
            stoichiometries.append(generate_stochiometry_dict(r['smiles']))
        except ValueError:
            continue


    S = np.zeros((len(molecules), len(reactions)))
    
    for i, s in enumerate(stoichiometries):

        if all([molecules_index.count(key) > 0 for key in s.keys()]):
            for key, value in s.items():
                
                S[molecules_index.index(key), i] = value
        else:
            continue
    
    S = pd.DataFrame(S, columns=list(map(lambda x: x['reaction_id'], reactions)), index=molecules)
    if not flag_too_big:
        st.write(S)
    else:
        st.download_button('Download!', data=S.to_csv(), file_name=upload_name.replace('.json', '.stoichiometry.csv'))


