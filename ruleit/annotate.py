from rdkit import Chem as chem 
from rdkit.Chem.Descriptors import ExactMolWt
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from rdkit.Chem import rdChemReactions as rdr
import requests
import urllib
import duckdb
import pandas as pd
import networkx as nx


class MoleculeBuffer:

    def __init__(self, file) -> None:
        """
        MoleculeBuffer allows to ease the annotation of chemical
        reaction networks by providing and storing annotations
        for molecules.
        """
        self.file = file

    def __transaction__(self, u):
        
        with duckdb.connect(self.file) as con:
            try:
                return con.sql(u).fetchall()
            except:
                return []
        # self.db = duckdb.connect(file)
    
    def create_table(self, overwrite=True):
        if overwrite:
            try:
                self.__transaction__("DROP TABLE molecules")
                
            except duckdb.duckdb.CatalogException:
                pass

        self.__transaction__(
            """
            CREATE TABLE molecules (
                key STRING, smiles STRING, inchi STRING, 
                inchikey STRING, title STRING, 
                cid INTEGER,
                mw FLOAT, formula STRING
            )
            """
        )
            


    def new_index(self):
        try:
            u = self.__transaction__("SELECT key FROM molecules")# .fetchall()
            n = max([int(item[0][1:]) for item in u]) + 1
            return 'm{:06d}'.format(n)
        except ValueError:
            return 'm000000'


    @staticmethod
    def compute_properties(molecule):

        u = dict(
            smiles=chem.MolToSmiles(molecule),
            inchi=chem.MolToInchi(molecule),
            inchikey=chem.MolToInchiKey(molecule)
        )
        try:
            molecule.UpdatePropertyCache()
        except:
            return u
        
        u['mw'] = ExactMolWt(molecule)
        u['formula'] = CalcMolFormula(molecule)
        return u
        
    def __qmol(self, inchikey):
        matches = self.__transaction__(
            """
            SELECT key, title, cid, inchikey,  smiles, formula, inchi, mw FROM molecules WHERE inchikey == '{:s}'
            """.format(inchikey)
        )
        return [dict(
            key=u[0], title=u[1], cid=u[2], inchikey=u[3], smiles=u[4], formula=u[5],
            inchi=u[6], mw=u[7]
        ) for u in matches]
        
    
    def query(self, q):
        """
        Queries must be specified using SMILES. MoleculeBuffer
        hands everything else
        """
        try:
            mol = chem.MolFromSmiles(q)
        except ValueError:
            raise IOError(f"unable to process {q}, check it is a SMILES")

        if mol is None:
            raise IOError(f"unable to process {q}, check it is a SMILES")
        
        inchikey = chem.MolToInchiKey(mol)
        u = self.__qmol(inchikey)
    
        if len(u) == 0:
            u = self.query_pubchem(q)
            if len(u) > 0:
                for match in u:
                    key = self.new_index()
                    self.__transaction__(
                        f"""
                        INSERT INTO molecules (key, title, cid, inchikey, smiles, formula, inchi, mw)
                            VALUES ('{key}', '{match['title']}', '{match['cid']}', '{match['inchikey']}', 
                            '{match['smiles']}', '{match['formula']}', '{match['inchi']}', '{match['mw']}')
                        """
                    )
                u = self.__qmol(inchikey)
            else:
                u = self.compute_properties(molecule=mol)
                key = self.new_index()
                self.__transaction__(
                    f"""
                    INSERT INTO molecules (key, inchikey, smiles, formula, inchi, mw)
                        VALUES ('{key}', '{u['inchikey']}', '{u['smiles']}', '{u['formula']}', '{u['inchi']}', '{u['mw']}')
                    """
                )
                u = self.__qmol(inchikey)
        
        return u


    @staticmethod
    def query_pubchem(smiles):
        """
        

        Parameters
        ----------
        smiles: str
            Smiles
        url: str
            URL. It must have a single placeholder to place the query.
        """
        url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/{:s}/property/CanonicalSMILES,Title,MolecularFormula,InChI,InChIKey,MolecularWeight/json?MaxRecords=5"
        url = urllib.parse.quote(url.format(smiles), safe=':/,?=')
        r = requests.get(url)
        out = []
        if r.status_code == 200:
            matches = r.json()['PropertyTable']['Properties']
            for mol in matches:
                try:
                    out.append(dict(
                        smiles=mol['CanonicalSMILES'],
                        inchi=mol['InChI'], inchikey=mol['InChIKey'], title=mol['Title'].lower(),
                        formula=mol['MolecularFormula'],
                        cid=mol['CID'], mw=mol['MolecularWeight'],
                    ))
                except KeyError:
                    print("unable to process {:s} match".format(smiles))
                    continue
        else:
            return []
        return out
 
    
def process_molecule_collection(collection, local_mol_db: MoleculeBuffer):
    molecule_dump = []
    stoich_dictionary = dict()
    unmatched_smiles = dict()
    unmatched_index = 0
    for mol in collection:
        
        smiles = chem.MolToSmiles(mol)
        try:
            m = local_mol_db.query(smiles)[0]
        except (IndexError, ValueError, IOError):
            try:
                key = unmatched_smiles[smiles]
            except KeyError:
                unmatched_smiles[smiles] = f'{unmatched_index:d}'
                unmatched_index += 1
                key = unmatched_smiles[smiles]

            m = dict(
                smiles=smiles,
                inchi="", inchikey="", title=smiles,
                formula="",
                cid="", mw="", key=key
            )
        
        molecule_dump.append(m)
        try:
            stoich_dictionary[m['key']] += 1
        except KeyError:
            stoich_dictionary[m['key']] = 1

    return molecule_dump, stoich_dictionary


def annotate_reaction_collection(reactions, local_mol_db):
    all_molecules = dict()
    links = dict()
    for entry in reactions:
        print(entry)
        reaction = rdr.ReactionFromSmarts(entry['smiles'], useSmiles=True)
        reaction_id = entry['reaction_id']
        reactants = reaction.GetReactants()
        products = reaction.GetProducts()
        reactant_stoich = dict()
        product_stoich = dict()

        reactants, reactant_stoich = process_molecule_collection(reactants, local_mol_db)
        products, product_stoich = process_molecule_collection(products, local_mol_db)
        

        for mol in reactants:
            
            all_molecules[mol['key']] = mol

        for mol in products:
            
            all_molecules[mol['key']] = mol

        
        for key, n in reactant_stoich.items():
            k = f"{reaction_id}:{key}"
            links[k] = dict(
                type="REACTS_IN", n=n
            )
            
            
        for key, n in product_stoich.items():
            k = f"{reaction_id}:{key}"
            links[k] = dict(
                type="REACTS_OUT", n=n
            )            

    return all_molecules, links


def generate_crn(reactions, molecules, links):

    G = nx.Graph()
    [r.update({'type': 'reaction'}) for r in reactions]
    [r.update({'label': r['reaction_id']}) for r in reactions]
    [r.update({'color': 'red'}) for r in reactions]

    [m.update({'type': 'molecule'}) for _, m in molecules.items()]
    [m.update({'label': m['title']}) for _, m in molecules.items()]
    [m.update({'color': 'blue'}) for _, m in molecules.items()]
    G.add_nodes_from([(r['reaction_id'], r) for r in reactions])
    G.add_nodes_from(list(molecules.items()))
    
    links = [(l.split(':')[0], l.split(':')[1], data) for l, data in links.items()]
    G.add_edges_from(links)
    return G