import streamlit as st
import pandas as pd
import json


st.title("Rule-It: Json 2 CSV 😳")


"""

We like to use JSON files to handle data internally, but it is not very handy for analysis. In this page,
you can convert the json file to set of csv files, which you can open later in Excel or Spreadsheets.

"""

upload = st.file_uploader(label="upload")
if upload:
    name = upload.name
    upload = json.load(upload)
    reactions = upload['reactions']
    molecules = upload['seeds']
    rules = upload['rules']

    
    
    

    
    buffer_reactions = pd.DataFrame.from_records(reactions)# .to_csv()
    buffer_molecules = pd.DataFrame.from_records(molecules)#.to_csv()
    buffer_rules = pd.DataFrame.from_records(rules)# .to_csv()

    cont1 = st.container(border=1)
    cont1.header('Reactions')
    col1, col2 = cont1.columns([6, 1.5])
    col1.metric('# reactions', len(buffer_reactions))
    col2.download_button(label='download!', key='download-reactions', data=buffer_reactions.to_csv(), file_name=name.replace('.json', 'reactions.csv'), type='primary')
    cont1.write(reactions[:5])

    cont2 = st.container(border=1)
    col1, col2 = cont2.columns([6, 1.5])
    col1.metric('# molecules', len(buffer_molecules))
    col2.download_button(label='download!', key='download-molecules', data=buffer_molecules.to_csv(), file_name=name.replace('.json', 'molecules.csv'), type='primary')
    cont2.write(molecules[:5])

    cont3 = st.container(border=1)
    col1, col2 = cont3.columns([6, 1.5])
    col1.metric('# rules', len(buffer_rules))
    col2.download_button(label='download!', key='download-rules', data=buffer_rules.to_csv(), file_name=name.replace('.json', 'rules.csv'), type='primary')
    cont3.write(rules[:5])
