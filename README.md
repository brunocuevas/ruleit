Welcome to Rule-It! This platform was created by Bruno Cuevas and Tyma Sokolskyi to enable researchers to study complex chemistry experiments, specifically, around the problem of abiogenesis and chemical evolution.

Our approach consists on:

1. Designing reaction rules.
2. Expanding a chemical reaction network using a collection of rules and a set of seed compounds.
3. Pruning the results using a linear-programming based algorithm to extract the minimal chemical reaction network.
4. Explore the results in an interactive and user friendly environment.

You can use this server online at https://ruleit.streamlit.app/

With this code you can also run Rule-it on your own machine, using Docker. Install Docker, clone this repository into your directory, and then run

docker build -t ruleit . 

to install Rule-it. To set up a local server run

docker run -p 8501:8501 ruleit

The preprint is at https://chemrxiv.org/engage/chemrxiv/article-details/671d306b1fb27ce124b0e6c9. Please cite as 
Cuevas-Zuviria, B. & Sokolskyi, T. (2024). Rule-it: an online platform for chemical reaction network explorations for chemical evolution. chemrxiv, doi: 10.26434/chemrxiv-2024-nqwpb
