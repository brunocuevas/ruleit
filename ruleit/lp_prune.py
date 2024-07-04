from rdkit import Chem as chem
from rdkit.Chem import rdChemReactions as rdr
from pulp import *

def split_reaction_string(reaction):
    reaction = rdr.ReactionFromSmarts(reaction['smiles'], useSmiles=True)
    #reactants, products = reaction['smiles'].split('>>')
    reactants = list(map(chem.MolToInchiKey, reaction.GetReactants()))
    products = list(map(chem.MolToInchiKey, reaction.GetProducts()))
    return reactants, products

def extract_molecules(reactions):
    molecules = []
    for reaction in reactions:
        mr, mp = reaction['smiles'].split('>>')
        mr = mr.split('.')
        mp = mp.split('.')
        for m in filter(lambda x: x not in molecules, mr + mp):
            molecules.append(m)
    return molecules


def process_reaction(reaction, reaction_name, 
                     molecule_reaction_map_reactant, molecule_reaction_map_product):

    
    reactants, products =  split_reaction_string(reaction)
    
    for r in filter(lambda x: x != "", reactants):
        try:
            molecule_reaction_map_reactant[r].append(reaction_name)
        except KeyError:
            molecule_reaction_map_reactant[r] = [reaction_name]

        

    for p in filter(lambda x: x != "", products):
        try:
            molecule_reaction_map_product[p].append(reaction_name)
        except KeyError:
            molecule_reaction_map_product[p] = [reaction_name]

    return None

    

def create_matrix(reactions, imports, exports):
    molecule_reaction_map_reactant = dict()
    molecule_reaction_map_product = dict()
    
    
    for i, reaction in enumerate(reactions):
        reaction_name = 'r{:06d}'.format(i)
        process_reaction(reaction, reaction_name, molecule_reaction_map_reactant, molecule_reaction_map_product)
    for i, reaction in enumerate(imports):
        reaction_name = 'i{:06d}'.format(i)
        process_reaction(reaction, reaction_name, molecule_reaction_map_reactant, molecule_reaction_map_product)
    for i, reaction in enumerate(exports):
        reaction_name = 'e{:06d}'.format(i)
        process_reaction(reaction, reaction_name, molecule_reaction_map_reactant, molecule_reaction_map_product)

    return molecule_reaction_map_reactant, molecule_reaction_map_product


def formalize_problem(reactions, seeds, hits):
    molecules = extract_molecules(reactions)
    import_reactions = [{"smiles": f">>{s}"} for s in seeds]
    export_reactions = [{"smiles": f"{h}>>"} for h in molecules]
    
    
    reaction_dict = dict()
    out_reactions = dict()
    hits = list(map(lambda x: chem.MolToInchiKey(chem.MolFromSmiles(x, sanitize=False)), hits))
    for i, r in enumerate(reactions):
        idx = 'r{:06d}'.format(i)
        reaction_dict[idx] = LpVariable(idx, 0, None, LpContinuous)
        out_reactions[idx] = r
    for i, s in enumerate(import_reactions):
        idx = 'i{:06d}'.format(i)
        reaction_dict[idx] = LpVariable(idx, 0, None, LpContinuous)
        out_reactions[idx] = s
    for i, (m, r) in enumerate(zip(molecules, export_reactions)):
        idx = 'e{:06d}'.format(i)
        
        m = chem.MolToInchiKey(chem.MolFromSmiles(m, sanitize=False))
        
        if m in hits:
            reaction_dict[idx] = LpVariable(idx, 0.1, None, LpContinuous)
            out_reactions[idx] = r
        else:
            reaction_dict[idx] = LpVariable(idx, 0.0, None, LpContinuous)
            out_reactions[idx] = r

    prob = LpProblem('MSP-Formose', LpMinimize)
    Mr, Mp = create_matrix(reactions, import_reactions, export_reactions)


    for m in molecules:
        m = chem.MolToInchiKey(chem.MolFromSmiles(m, sanitize=False))
        try:
            product_of = Mr[m]
        except KeyError:
            product_of = []
        try:
            reactant_of = Mp[m]
        except KeyError:
            reactant_of = []
        prob += lpSum(
            [(1.0 * reaction_dict[x]) for x in product_of] + 
            [(-1.0 * reaction_dict[x]) for x in reactant_of]
        ) == 0.0

    prob += lpSum([(1.0 * reaction) for key, reaction in reaction_dict.items() if key[0] == 'r'])
    prob.solve()
    for r in prob.variables():
        out_reactions[r.name]['active'] = (r.varValue != 0.0)
    return prob.sol_status, out_reactions

# for v in prob.variables():
#     print(v.name, v.varValue)