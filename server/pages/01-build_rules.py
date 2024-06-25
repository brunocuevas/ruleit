import streamlit as st
from rdkit.Chem import rdChemReactions as rdr
from rdkit.Chem import Draw
from rdkit import Chem as chem

st.set_page_config(
    page_title="Rule-it:Build",
    page_icon="",
)


st.title('Rule-it: Reaction rule builder 🏗️')

st.markdown(
    """
    Making reaction rules can be quite painful. In this app,
    we seek to simplify it. In the first box, you can include
    a set of compounds that you want to see reaction (or not!). In
    the second box, you write your SMARTS, so you can visualize the final
    result.
    """
)

reactants = st.text_area('Introduce your reactants here')
reaction_rule = st.text_area('Introduce your expression here')



reactants = [u.split('.') for u in reactants.split('\n')]

if len(reactants) > 0 and len(reaction_rule) > 0:
    st.header('Results')
    chemical_reaction = rdr.ReactionFromSmarts(reaction_rule)
    n_reactants = chemical_reaction.GetNumReactantTemplates()
    flag = False
    for r in reactants:
        r = [chem.MolFromSmiles(x) for x in r]
        try:
            products = chemical_reaction.RunReactants(r)
        except ValueError:
            continue
        if len(products) > 0:
            for p in products:
                flag = True
                st.image(Draw.ReactionToImage(
                    rdr.ReactionFromSmarts(
                        '>>'.join(
                            [
                                '.'.join([chem.MolToSmiles(m) for m in r]),
                                '.'.join([chem.MolToSmiles(m) for m in p])
                            ]
                        ), useSmiles=True
                    )
                ))
    if not flag:
        st.write('Not found!')

