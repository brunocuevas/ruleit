from rdkit import Chem as chem
from rdkit.Chem import rdChemReactions as rdr
import hashlib
import base64


def hash_molecules(smiles: str = None, mol: str = None) -> str:
    # mol = None
    if mol is None:
        mol = chem.MolFromSmiles(smiles)
    elif smiles is None:
        mol = chem.Mol(mol)
    if mol is None:
        return None
        
    chem.RemoveStereochemistry(mol)
    inchikey = chem.MolToInchiKey(mol)
    if inchikey == "":
        raise RuntimeError("couldn't generated identifier - {:s}".format(chem.MolToSmiles(mol)))
    if mol.GetNumHeavyAtoms() > 6:
        return 'mol-l-{:s}'.format(inchikey[:-2])
    else:
        return 'mol-s-{:s}'.format(inchikey)


def hash_reactions(smiles: str) -> str:
    """
    Generates a hash out of a chemical reaction in order to recognize
    their identity regardless of the SMILES employed

    """
    reaction = rdr.ReactionFromSmarts(smiles, useSmiles=True)
    reactants = reaction.GetReactants()
    products = reaction.GetProducts()
    reactant_identifiers = list(map(lambda x: hash_molecules(mol=x), reactants))
    product_identifiers = list(map(lambda x: hash_molecules(mol=x), products))

    merged_ids = ','.join(sorted(reactant_identifiers)) + '->' + ','.join(sorted(product_identifiers))
    hasher = hashlib.sha1(merged_ids.encode('utf-8')).digest()[:12]
    new_key = 'rxn-' + base64.urlsafe_b64encode(hasher).decode('utf-8')
    return new_key