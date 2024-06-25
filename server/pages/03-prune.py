import streamlit as st
import json
from ruleit.lp_prune import formalize_problem
from pulp import *
import pandas as pd
from io import StringIO, BytesIO
import networkx as nx 


st.title("Rule-it: Network prunning ✂️")

st.markdown(
    """
    Chemical reaction networks can be extremely complex
    and cumbersome. Reactions will take place depending on
    parameters that we might be unable to control or to know
    at all, giving place to infinity of different results.


    We seek to study these complex chemical reaction networks
    by prunning the results of the reaction-rule expansion by
    extracting the minimal reaction network that leads to the
    observed results. We do so by using linear-programming.
    """
)

st.markdown("Introduce the file that you obtained from the previous step")
u = st.file_uploader(label="network")
st.markdown("Introduce the list of smiles that you used in the previous step")
v = st.text_area(label="seeds")
st.markdown("Introduce a SMILES list of the detected compounds")
w = st.text_area(label="hits")

if u and v:

    u = json.load(u)
    v = [line.strip() for line in v.split('\n')]
    w = [line.strip() for line in w.split('\n')]

    status, results = formalize_problem(
        reactions=u['reactions'],
        seeds=v, 
        hits=w
    )
    if status == 1:
        results = dict((key, item) for key, item in filter(lambda x: x[0][0] == 'r', results.items()))
        results = pd.DataFrame.from_records(results).T
        st.metric(label="Coverage", value=len(results.query('active == True'))/len(results))
        st.dataframe(results.query('active == True'), use_container_width=True)
        pruned_reactions = results.query('active == True')
        # st.write(pruned_reactions.sample(n=min([25, len(pruned_reactions)])))
        u['reactions'] = pruned_reactions.to_dict(orient='records')
        st.download_button(
            label="Download results", 
            data=json.dumps(u, indent=4), 
            file_name="pruned-expansion.json"
        )

    else:
        st.warning("unable to solve LP problem!")
# for v in prob.variables():
#     print(v.name, v.varValue)
    

