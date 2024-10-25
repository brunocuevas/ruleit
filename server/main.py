import streamlit as st


st.set_page_config(
    page_title="Rule-it:Wellcome",
    page_icon="👋",
)


st.write("# Welcome to Rule-it!👋")
st.sidebar.success("Select a demo above.")

st.markdown(
    """
    Wellcome to Rule-It! This platform was created by Bruno Cuevas
    and Tyma Sokolskyi to enable researchers to study complex
    chemistry experiment, specifically, around the problem of
    abiogenesis and chemical evolution.


    Our approach consists on:
    1. Creating reaction rules.
    2. Expanding a chemical reaction network using a collection of rules 
    and a set of seed compounds.
    3. Prune the results using a linear-programming based algorithm to extract
    the minimal chemical reaction network.
    4. Explore the results in an interactive and user friendly environment.


    You can also run this server on your own machine. Check the 
    code at https://github.com/brunocuevas/ruleit

    """
)
