import streamlit as st
import json
from ruleit.annotate import MoleculeBuffer
import pandas as pd
from stqdm import stqdm

if 'direction' not in st.session_state:
    st.session_state['direction'] = True


@st.cache_resource()
def local_buffer():
    buffer = MoleculeBuffer("moleculebuffer.db")
    buffer.create_table(overwrite=False)
    return buffer

st.title("Rule-It: Map SMILES 🗺️")

"""
Introduce a list of names here, and we will try to map them
to their smiles through the Pubchem API.
"""

if st.toggle("Invert: Smiles to name"):
    st.session_state['direction'] = False

if st.session_state["direction"]:
    names = st.text_area(label="introduce your {:s} here".format('names'))
else:
    names = st.text_area(label="introduce your {:s} here".format('smiles'))

process_button = st.button(label='run')

if process_button and names:

    hits = []

    buffer = local_buffer()

    names = list(map(lambda x: x.strip(), names.split('\n')))

    for name in stqdm(names):
        if st.session_state['direction']:
            try:
                u = buffer.query_pubchem(q=name, use='name')[0]
            except IndexError:
                continue
        else:
            try:
                u = buffer.query_pubchem(q=name, use='smiles')[0]
            except IndexError:
                continue

        u['query'] = name

        hits.append(u)

    hits = pd.DataFrame.from_records(hits)

    col1, col2 = st.columns(2)

    col1.metric(
        label="# succesful hits", value=len(hits.dropna(subset=['smiles']))
    )
    col2.metric(
        label="succes rate", value=100 * len(hits.dropna(subset=['smiles'])) / len(names)
    )

    st.dataframe(
        hits[['query', 'title', 'smiles', 'mw', 'formula']],
        use_container_width=True
    )

