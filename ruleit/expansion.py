import yaml
from rdkit.Chem import rdChemReactions as rdr
from rdkit.Chem import AllChem as achem
from rdkit import Chem as chem
from rdkit.Chem.Descriptors import MolWt
import itertools
import numpy as np
from random import shuffle



def cartessian_product(seed_set, n):
    u = [seed_set] * n
    return itertools.product(*u)


def generate_smiles(s, p):
    return '>>'.join(
        ['.'.join([chem.MolToSmiles(m) for m in s]),
        '.'.join([chem.MolToSmiles(m) for m in p])]
    )


def _expansion(seeds, reaction_rules, max_reactions=1000):
    """

    It performs a simple single iteration update of the
    pool of molecules. It returns the discovered reactions, their
    classification, and a pool of new molecules.

    """

    
    seeds = list(map(lambda x: x.strip(), seeds))
    seeds = list(set(seeds))
    
    seed_molecules = []
    output = dict()
    
    for m in seeds:
        try:
            mol = chem.MolFromSmiles(m)
            if mol is None:
                
                continue
            seed_molecules.append(mol)
        except ValueError:
            
            pass
    

    discovered_reactions = []
    discovered_substrates = []
    discovered_substrates = discovered_substrates + seeds
    idx = 0
    for reaction in reaction_rules['reactions']:

        chemical_reaction = rdr.ReactionFromSmarts(reaction['smarts'])
        n_reactants = chemical_reaction.GetNumReactantTemplates()
        chemical_reaction.Initialize()
        subset_reactants = [m for m in seed_molecules if chemical_reaction.IsMoleculeReactant(m)]
        substrates = cartessian_product(subset_reactants, n_reactants)
        
        for s in substrates:
            p = chemical_reaction.RunReactants(s)
            if len(p) > 0:
                for batch in p:
                    discovered_reactions.append(
                        dict(
                            smiles=generate_smiles(s, batch),
                            rule=reaction['name'], 
                            reaction_id=f"r{idx:06d}"
                        )                        
                    )
                    idx += 1
                    for m in batch:
                        discovered_substrates.append(
                            chem.MolToSmiles(m)
                        )
                    if idx >= max_reactions:
                        print("max reactions reached")
                        output['discovered-molecules'] = _prune_molecules(discovered_substrates)
                        output['discovered-reactions'] = discovered_reactions
                        return output
                    
    
    output['discovered-molecules'] = _prune_molecules(discovered_substrates)
    output['discovered-reactions'] = discovered_reactions

    return output


def _prune(discovered_reactions):
    out = []
    index = []
    for r in discovered_reactions:
        if r['smiles'] not in index:
            out.append(r)
            index.append(r['smiles'])
    return out

def _prune_molecules(discovered_molecules):
    out = []
    for mol in discovered_molecules:
        if mol not in out:
            out.append(mol)
    out.sort(key=len)
    return out

def probablistic_expansion(seeds, reaction_rules, rule_probability, iterations=1000, conditions=None):
    """

    Performs a probablistic expansion, where the rule_probability parameters
    determine how likely is to apply one rule or the other. Everytime a reaction
    is applied, a random subset of its reactants is generated, and then only one reaction
    is carried out.

    """

    reaction_rules_copy = dict(
        reactions=[r.copy() for r in reaction_rules['reactions']]
    )
    seeds = list(map(lambda x: x.strip(), seeds))
    seeds = list(set(seeds))
    
    seed_molecules = []
    output = dict()
    
    for m in seeds:
        try:
            mol = chem.MolFromSmiles(m)
            if mol is None:
                
                continue
            seed_molecules.append(mol)
        except ValueError:
            
            pass
    

    discovered_reactions = []
    discovered_substrates = []
    discovered_substrates = discovered_substrates + seeds
    idx = 0

    for reaction in reaction_rules_copy['reactions']:
        x = reaction['smarts']
        chemical_reaction = rdr.ReactionFromSmarts(x)
        chemical_reaction.Initialize()
        reaction['reaction_object'] = chemical_reaction
        reaction['number_reactants'] = chemical_reaction.GetNumReactantTemplates()
        reaction['reactants'] = dict()
        flag = True
        for i in range(reaction['number_reactants']):
            template = chemical_reaction.GetReactantTemplate(i)
            reaction['reactants'][f'role_{i:02d}'] = [m for m in seed_molecules if m.GetSubstructMatch(template)]
            if len(reaction['reactants'][f'role_{i:02d}']) == 0:
                flag = False
        reaction['active'] = flag

    for i in range(iterations):
    
        reaction = np.random.choice(reaction_rules_copy['reactions'], p=rule_probability, size=1)[0]
        chemical_reaction = reaction['reaction_object']

        if not reaction['active']:
            continue
        
        s = []

        for _, roles in reaction['reactants'].items():
            m = np.random.choice(roles, p=np.ones(len(roles)) / len(roles), size=1)[0]
            s.append(m)
        
        p = chemical_reaction.RunReactants(s)

        if len(p) > 0:
            for batch in p:
                if conditions is not None:
                    if not safety_check(batch, conditions=conditions):
                        continue

                discovered_reactions.append(
                    dict(
                        smiles=generate_smiles(s, batch),
                        rule=reaction['name'], 
                        reaction_id=f"r{idx:06d}"
                    )                        
                )
                idx += 1
                for m in batch:
                    discovered_substrates.append(
                        chem.MolToSmiles(m)
                    )



    output['discovered-molecules'] = _prune_molecules(discovered_substrates)
    output['discovered-reactions'] = discovered_reactions

    

    return output


def iterative_probabilistic_expansion(seeds, reaction_rules, iterations, rounds, conditions=None):
    expansion_metrics = []
    for i in range(rounds):
        current_expansion = probablistic_expansion(
            seeds, reaction_rules=dict(reactions=reaction_rules),
            rule_probability=np.ones(len(reaction_rules)) / len(reaction_rules), 
            iterations=iterations,
            conditions=conditions
        )
        
        try:
            out['discovered-reactions'] = _prune(out['discovered-reactions'] + current_expansion['discovered-reactions'])
            out['discovered-molecules'] = current_expansion['discovered-molecules']
            
        except NameError:
            out = current_expansion
            out['discovered-reactions'] = _prune(out['discovered-reactions'])

        seeds = out['discovered-molecules']

        nr = len(out['discovered-reactions'])
        nm = len(out['discovered-molecules'])

        expansion_metrics.append(
            dict(iteration=i, n_reactions = nr, n_molecules = nm)
        )

    for i, reaction in enumerate(out['discovered-reactions']):
        reaction['reaction_id'] = f'r{i:06d}'
    
    output_content = dict(
        rules=reaction_rules, 
        seeds=seeds,
        reactions=out['discovered-reactions'],
        metrics=expansion_metrics
    )
    return output_content


def safety_check(products, conditions):
    for product in products:
        product = chem.MolFromInchi(chem.MolToInchi(product))
        if product is None:
            return False
        cv = check_valence(product, conditions)
        cm = check_mass(product, conditions)
        if cv and cm:
            continue
        else:
            return False
    return True

def check_valence(molecule, conditions):
    # chem.SanitizeMol(molecule)
    
    for atom in molecule.GetAtoms():

        try:
            valence = atom.GetTotalValence()
        except RuntimeError:
            return False
        
        symbol = atom.GetSymbol()
        
        try:
            allowed_valence = conditions['valence'][symbol] 
        except KeyError:
            continue

        if not (valence in allowed_valence):
            return False
        
    return True

def check_mass(molecule, conditions):
    try:
        mw = MolWt(molecule)
    except RuntimeError:
        return False
    if mw > conditions['mass']:
        return False
    else:
        return True