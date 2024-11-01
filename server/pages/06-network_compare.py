import streamlit as st
from ruleit.hash import hash_molecules, hash_reactions
from ruleit.annotate import annotate_reaction_collection, MoleculeBuffer, generate_crn
import pandas as pd
import networkx as nx
import json


@st.cache_data()
def upload_reaction_network(reactions, _buffer):
    molecules, links = annotate_reaction_collection(reactions=reactions, local_mol_db=_buffer)
    return dict(reactions=reactions, molecules=molecules, links=links)

@st.cache_resource()
def local_buffer():
    buffer = MoleculeBuffer("moleculebuffer.db")
    buffer.create_table(overwrite=False)
    return buffer


def compute_molecule_overlap(seeds1, seeds2):
    
    seeds1 = [hash_molecules(smiles=mol['smiles']) for mol in seeds1 if mol['smiles'] != ""]
    seeds2 = [hash_molecules(smiles=mol['smiles']) for mol in seeds2 if mol['smiles'] != ""]
    overlapping_seeds_1 = [s for s in seeds1 if s in seeds2]
    overlapping_seeds_2 = [s for s in seeds2 if s in seeds1]
    n_overlap = min(len(overlapping_seeds_1), len(overlapping_seeds_2))
    return n_overlap, n_overlap / len(list(set(seeds1 + seeds2)))


def compute_reaction_overlap(reactions1, reactions2):
    reactions1 = [hash_reactions(reaction['smiles']) for reaction in reactions1]
    reactions2 = [hash_reactions(reaction['smiles']) for reaction in reactions2]
    overlapping_reactions_1 = [s for s in reactions1 if s in reactions2]
    overlapping_reactions_2 = [s for s in reactions2 if s in reactions1]
    n_overlap = min([len(overlapping_reactions_1), len(overlapping_reactions_2)])
    return n_overlap, n_overlap / len(list(set(reactions1 + reactions2)))



def get_graph_attributes(G):
    btw = nx.betweenness_centrality(G, normalized=True)
    eig = nx.eigenvector_centrality(G)
    prnk = nx.pagerank(G)
    
    for node, node_data in filter(lambda x: x[1]['type'] == 'molecule', G.nodes(data=True)):

        node_data['btw'] = btw[node]
        node_data['eig'] = eig[node]
        node_data['prnk'] = prnk[node]
        node_data['deg'] = G.degree[node]
        node_data['hash'] = hash_molecules(smiles=node_data['smiles'])

    for node, node_data in filter(lambda x: x[1]['type'] == 'reaction', G.nodes(data=True)):

        node_data['btw'] = btw[node]
        node_data['eig'] = eig[node]
        node_data['prnk'] = prnk[node]
        node_data['deg'] = G.degree[node]
        node_data['hash'] = hash_reactions(node_data['smiles'])

    return G


def extract_molecule_attributes(G):
    return pd.DataFrame.from_records(
        list(map(lambda x: x[1], filter(lambda x: x[1]['type'] == 'molecule', G.nodes(data=True))))
    )[['hash', 'title', 'formula', 'mw', 'deg', 'prnk', 'eig', 'btw', 'smiles']]

def extract_reaction_attributes(G):
    return pd.DataFrame.from_records(
        list(map(lambda x: x[1], filter(lambda x: x[1]['type'] == 'reaction', G.nodes(data=True))))
    )[['hash',  'deg', 'prnk', 'eig', 'btw', 'smiles']]



st.title("Rule-It: Compare 😳")

"""
This tool enables you to compare the output of two different runs (e.g. expansions with two different
seed sets).
"""

data = False
col1, col2 = st.columns(2)
upload_1 = col1.file_uploader(label="upload-1")
upload_2 = col2.file_uploader(label="upload-2")

if upload_1 and upload_2:
    buffer = local_buffer()
    data = dict(
        upload_1=upload_reaction_network(reactions=json.load(upload_1)['reactions'], _buffer=buffer),
        upload_2=upload_reaction_network(reactions=json.load(upload_2)['reactions'], _buffer=buffer),
    )
    G1 = generate_crn(**data['upload_1'])
    G2 = generate_crn(**data['upload_2'])
    cont = st.container(border=1)

    cont.header("Molecules")

    n_overlapping_molecules, p_overlapping_molecules = compute_molecule_overlap(
        seeds1=list(data['upload_1']['molecules'].values()), 
        seeds2=list(data['upload_2']['molecules'].values())
    )

    col1, col2 = cont.columns(2)
    col1.metric('# overlapping hits', value=n_overlapping_molecules)
    col2.metric(r'% overlapping hits', value=p_overlapping_molecules)

    G1 = get_graph_attributes(G1)
    G2 = get_graph_attributes(G2)

    mol_attributes_1 = extract_molecule_attributes(G1)
    mol_attributes_2 = extract_molecule_attributes(G2)
    

    u = pd.merge(mol_attributes_1, mol_attributes_2, on=['hash', 'title', 'formula'], how='inner')

    mol_attr_1 = st.selectbox('show_attr', options=['deg', 'prnk', 'eig', 'btw'])

    st.scatter_chart(
        data=u, x=mol_attr_1 + '_x', y=mol_attr_1 + '_y', color='title'
    )


    cont = st.container(border=1)

    cont.header("Reactions")
    n_overlapping_reactions, p_overlapping_reactions = compute_reaction_overlap(
        reactions1=data['upload_1']['reactions'], reactions2=data['upload_2']['reactions']
    )

    col1, col2 = cont.columns(2)
    col1.metric('# overlapping hits', value=n_overlapping_reactions)
    col2.metric(r'% overlapping hits', value=p_overlapping_reactions)


    reaction_attributes_1 = extract_reaction_attributes(G1)
    reaction_attributes_2 = extract_reaction_attributes(G2)

    u = pd.merge(reaction_attributes_1, reaction_attributes_2, on=['hash'], how='inner')

    mol_attr_2 = st.selectbox('show_attr-2', options=['deg', 'prnk', 'eig', 'btw'])

    st.scatter_chart(
        data=u, x=mol_attr_2 + '_x', y=mol_attr_2 + '_y', color='smiles_x'
    )



    

