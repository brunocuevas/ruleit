import streamlit as st
import json
from ruleit.annotate import annotate_reaction_collection, MoleculeBuffer, generate_crn
import numpy as np
import networkx as nx
import plotly.express as px
import pandas as pd
from streamlit_agraph import agraph, Node, Edge, Config



@st.cache_data()
def upload_reaction_network(reactions, _buffer):
    molecules, links = annotate_reaction_collection(reactions=reactions, local_mol_db=_buffer)
    return dict(reactions=reactions, molecules=molecules, links=links)

@st.cache_resource()
def local_buffer():
    buffer = MoleculeBuffer("moleculebuffer.db")
    buffer.create_table(overwrite=False)
    return buffer

st.title("Rule-It: Explore 😳")

"""
This page enables you to enrich the output of the previous runs. Be aware that it requires making requests
to other databases, so it will take a while (or forever) to annotate very big reaction networks.
"""

data = False
upload = st.file_uploader(label="upload")
if upload:
    buffer = local_buffer()
    data = upload_reaction_network(reactions=json.load(upload)['reactions'], _buffer=buffer)

if data:
    

    G = generate_crn(**data)

    st.subheader("Network size")

    col1, col2 = st.columns(2)
    col1.metric("Number of nodes:", G.number_of_nodes())
    col2.metric("Number of edges:", G.number_of_edges())



    st.subheader("Degree")

    
    btw = nx.betweenness_centrality(G, normalized=True)
    eig = nx.eigenvector_centrality(G)
    prnk = nx.pagerank(G)
    
    for node, node_data in filter(lambda x: x[1]['type'] == 'molecule', G.nodes(data=True)):

        node_data['btw'] = btw[node]
        node_data['eig'] = eig[node]
        node_data['prnk'] = prnk[node]
        node_data['deg'] = G.degree[node]
    
    
    molecular_data = pd.DataFrame.from_records(
        list(map(lambda x: x[1], filter(lambda x: x[1]['type'] == 'molecule', G.nodes(data=True))))
    )[['key', 'title', 'formula', 'mw', 'deg', 'prnk', 'eig', 'btw', 'smiles']]

    

    col1, col2 = st.columns(2)
    col1.metric("Average degree", np.mean(molecular_data['deg']))
    col2.metric("Average density", np.round(nx.density(G), 3))
    
    molecule_degree = [item[1] for item in list(G.degree)]
    
    fig = px.histogram(data_frame=molecular_data, x='deg', labels={"x": "Molecule Degree"})
    event = st.plotly_chart(fig, key="deg", use_container_width=True)
    st.write(molecular_data.sort_values(by='deg', ascending=False))

    st.caption(
        """
        **mw**: molecular weight; 
        **deg**: degree; 
        **prnk**: page-rank; 
        **eig**: eigenvector centrality; 
        **btw**: betweenness centrality
        """
    )

    nodes = [Node(id=node, **node_data) for node, node_data in G.nodes(data=True)]
    edges = [Edge(source=e1, target=e2, **edge_data) for e1, e2, edge_data in G.edges(data=True)]
    enable_physics = len(nodes) < 200
    config = Config(width=600,
                height=600,
                directed=False, 
                physics=enable_physics, 
                hierarchical=False,
                # **kwargs
                )

    return_value = agraph(nodes=nodes, 
                        edges=edges, 
                        config=config)