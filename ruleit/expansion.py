import yaml
from rdkit.Chem import rdChemReactions as rdr
from rdkit import Chem as chem
import itertools


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
        substrates = cartessian_product(seed_molecules, n_reactants)
        
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
                        output['discovered-molecules'] = discovered_substrates
                        output['discovered-reactions'] = discovered_reactions
                        return output
                    
    
    output['discovered-molecules'] = discovered_substrates
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