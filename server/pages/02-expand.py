from ruleit.expansion import iterative_probabilistic_expansion
import streamlit as st 
from rdkit.Chem import rdChemReactions as rdr
from rdkit.Chem import Draw
from rdkit import Chem as chem
import networkx as nx
import json
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import random



if 'reaction-rules' not in st.session_state:
    st.session_state['reaction-rules'] = ""
if 'seeds' not in st.session_state:
    st.session_state['seeds'] = ""


st.title("Rule-it: Network expansion 🥧")
st.header("Summary")
"""
RuleIt is a platform aimed to the generation of chemical reaction networks
using three sources of knowledge: i) chemical reaction rules; 2) databases
of chemical reactions; and iii) artificial intelligence.
"""

"""

The network generator allows you to obtain a chemical reaction network
from a set of rules and seed compounds. For that, we use a very simple algorithm:
the network expansion.

"""

"""
First, using your starting seed set, we compute every chemical reaction
possible using that set. Then, we compute the products of those reactions. Finally, 
we repeat until the maximum number of iterations (right now, only 3).

"""

st.header("Let's rule it")

"""

First, we need the maximum number of iterations that you are going to perform.
As this number grows very fast, we are limiting it to three iterations. 

"""

col1, col2, col3 = st.columns(3)
n_expansions = col1.slider('# expansion iterations', value=4, min_value=1, max_value=5)
iterations_per_expansion = col2.slider("# number reactions per expansion", min_value=1000, max_value=50000, value=10000, step=1000)
maximum_mass = col3.slider('maximum mass', value=500, min_value=250, max_value=1000)
seed = st.text_input('random-seed', value=42)
np.random.seed(int(seed))

"""
Next, we need the reaction rules. Reaction rules can become quite 'mind-bending', we are aware. That is why we
provide a reaction rule builder so you can debug your own reaction rules. We recomend you that you simply start
some spreadsheet with your reaction rules, so you can copy them here as plain text. Rules will be number by the
order that they appear. We use Python, so we count from zero (your first rule, is our rule r000).
"""
example_rules = st.button("give me an example",  key="example-1")

if example_rules:
    st.session_state['reaction-rules'] = "[CD3:1](=[OD1:2])[O:3].[ND1:4]-[CD1,CD2:5]>>[C:1](=[O:2])[N:4][C:5].[O:3]"
    st.session_state['seeds'] = "CCC(=O)O\nNCC(=O)O"
    

reaction_rules = st.text_area(
    label="reaction rules", key='reaction-rules'
)

"""
Now, we need the seeds. Seeds are the compounds that you use to start your network. In a complex chemical
setting, those would be the reactants that you introduce, but you might also want to consider these as the
starting points for the network that you might generate here.
"""


seeds = st.text_area(
    "seeds", key="seeds"
)

"""
By clicking process, you will start the creation of the reaction network. Notice that, depending on the number
of iterations, the size of the seed and reaction rule set, and the intrinsic nature of the reaction rules that you
introduce, this step can be either instantaneously or take forever. 

"""

conditions = dict(
    mass=maximum_mass,
    valence=dict(
        C=[4],
        N=[3],
        O=[2]
    )
)


if st.button('process'):
    
    st.header("Results")
    st.divider()

    reaction_rules = reaction_rules.split('\n')
    seeds = seeds.split('\n')
    reaction_rules = [dict(name='r{:06d}'.format(i), smarts=r) for i, r in enumerate(reaction_rules)]

    outfile = iterative_probabilistic_expansion(
        seeds=seeds, 
        reaction_rules=reaction_rules, 
        iterations=iterations_per_expansion, 
        rounds=n_expansions, 
        conditions=conditions
    )

    reactions = pd.DataFrame.from_records(outfile['reactions'])
    reaction_counts = reactions.groupby('rule', as_index=False).count().rename(columns={'smiles':'#number-reactions'})

    
    st.download_button(label='Download', data=json.dumps(outfile, indent=4), file_name='expansion.json', type="primary")

    st.bar_chart(data=reaction_counts, x='rule', y='#number-reactions')

    col1, col2 = st.columns(2)
    col1.bar_chart(
        data=pd.DataFrame.from_records(outfile['metrics']), x='iteration', y="n_reactions"
    )
    col2.bar_chart(
        data=pd.DataFrame.from_records(outfile['metrics']), x='iteration', y="n_molecules"
    )
    col1, col2 = st.columns(2)
    col1.metric("Number of compounds", len(outfile['seeds']))
    col2.metric("Number of reactions", len(outfile['reactions']))

    
    st.divider()
    """
    Some of the reactions found in the expansion
    """    
    display_reactions = outfile['reactions'].copy()
    random.shuffle(display_reactions)
    display_reactions = display_reactions[:10]

    st.divider()


    st.subheader("Reaction examples")
    for r in display_reactions:
        st.image(Draw.ReactionToImage(
                rdr.ReactionFromSmarts(
                    r['smiles'], useSmiles=True
                )
            ))
        st.write(r)
