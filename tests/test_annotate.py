import unittest
from ruleit.annotate import MoleculeBuffer, annotate_reaction_collection, generate_crn
import json
import random


class TestAnnotate(unittest.TestCase):

    def test_buffer(self):
        mb = MoleculeBuffer(file="tests/mol-buffer.db")
        
        mb.create_table()
        
        mb.query("CCO")
        u = mb.query("CCO")
        u = mb.query("CCO")
        self.assertEqual(len(u), 1)
        self.assertEqual(u[0]['key'], 'm000000')


        u = mb.query("NCCO")
        self.assertEqual(len(u), 1)
        self.assertEqual(u[0]['key'], 'm000001')

    def test_annotate_collection(self):
        mb = MoleculeBuffer(file="tests/mol-buffer.db")
        mb.create_table()
        u = [{
            "reaction_id": "r000001",
            "smiles": "CCO>>CC=O"
        }, {
            "reaction_id": "r000001",
            "smiles": "CC=O.O>>CC(=O)O"
        }]
        molecules, links = annotate_reaction_collection(
            reactions=u, local_mol_db=mb
        )

        alcohol = mb.query("CCO")
        aldehyde = mb.query("CC=O")
        carboxyl = mb.query("CC(=O)O")

        self.assertEqual(alcohol[0]['key'], 'm000000')
        self.assertEqual(aldehyde[0]['key'], 'm000001')
        self.assertEqual(carboxyl[0]['key'], 'm000003')


        self.assertEqual(len(molecules.keys()), 4)

    def test_exception_1(self):
        mb = MoleculeBuffer(file="tests/mol-buffer.db")
        mb.create_table()
        annotate_reaction_collection(
            [{
                "smiles": "NCC(=O)NCC(=O)O.NCC(=O)NCC(=O)O>>NCC(=O)NCC(=O)NCC(=O)NCC(=O)O.O",
                "rule": "r000000",
                "reaction_id": "r000000"
            },
            {
                "smiles": "NCC(=O)NCC(=O)O.NCC(=O)NCC(=O)NCC(=O)NCC(=O)O>>NCC(=O)NCC(=O)NCC(=O)NCC(=O)NCC(=O)NCC(=O)O.O",
                "rule": "r000000",
                "reaction_id": "r000001"
            }], local_mol_db=mb
        )
        m = mb.query('NCC(=O)NCC(=O)O')
        self.assertEqual(len(m), 1)
        
        
    def test_network_generation(self):
        mb = MoleculeBuffer(file="tests/mol-buffer.db")
        mb.create_table()
        reactions = [{
            "reaction_id": "r000001",
            "smiles": "CCO>>CC=O"
        }, {
            "reaction_id": "r000002",
            "smiles": "CC=O.O>>CC(=O)O"
        }]
        molecules, links = annotate_reaction_collection(reactions, mb)
        G = generate_crn(reactions, molecules, links)
        self.assertEqual(G.number_of_nodes(), 6)

    def test_large_network(self):
        mb = MoleculeBuffer(file="tests/mol-buffer.db")
        mb.create_table(overwrite=True)
        with open('tests/expansion-2it.json', 'r') as f:
            reactions = json.load(f)['reactions']
        # random.shuffle(reactions)
        molecules, links = annotate_reaction_collection(reactions[:50], mb)
        
    def test_large_network_2(self):
        mb = MoleculeBuffer(file="tests/mol-buffer.db")
        mb.create_table(overwrite=True)
        with open('tests/expansion-1.json', 'r') as f:
            reactions = json.load(f)['reactions']
        # random.shuffle(reactions)
        molecules, links = annotate_reaction_collection(reactions[:50], mb)

    def test_large_network_3(self):
        mb = MoleculeBuffer(file="tests/mol-buffer.db")
        mb.create_table(overwrite=True)
        with open('tests/pruned-expansion-all.json', 'r') as f:
            reactions = json.load(f)['reactions']
        # random.shuffle(reactions)
        molecules, links = annotate_reaction_collection(reactions, mb)
        crn = generate_crn(reactions, molecules, links)


if __name__ == "__main__":
    unittest.main()

