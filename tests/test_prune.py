from ruleit.lp_prune import *
from ruleit.expansion import iterative_probabilistic_expansion
import unittest


class TestLP(unittest.TestCase):

    def test_lp_prune_1(self):

        status, results = lp_prune(
            reactions=[
                {"smiles":"C>>CC", "reaction_id":"r000000"},
                {"smiles":"CC>>CCC", "reaction_id":"r000001"},
                {"smiles":"CCC>>CCN", "reaction_id":"r000002"},
                {"smiles":"CCC>>CCCC", "reaction_id":"r000003"}
            ], seeds=["C"],
            hits=["CCN"]
        )
        self.assertEqual(results['r000000']['active'], True)
        self.assertEqual(results['r000001']['active'], True)
        self.assertEqual(results['r000002']['active'], True)
        self.assertEqual(results['r000003']['active'], False)

    def test_lp_prune_2(self):

        status, results = lp_prune(
            reactions=[
                {"smiles":"C(O)=O>>CC", "reaction_id":"r000000"},
                {"smiles":"CC>>CCC", "reaction_id":"r000001"},
                {"smiles":"CCC>>CCN", "reaction_id":"r000002"},
                {"smiles":"CCC>>CCCC", "reaction_id":"r000003"}

            ], seeds=["OC=O"],
            hits=["NCC", "CCCC"]
        )
        if status == 1:
            self.assertEqual(results['r000000']['active'], True)
            self.assertEqual(results['r000001']['active'], True)
            self.assertEqual(results['r000002']['active'], True)
            self.assertEqual(results['r000003']['active'], True)
            
        else:
            self.assertEqual(status, -1)

    def test_lp_prune_3(self):

        status, results = lp_prune(
            reactions=[
                {"smiles": "CC=O.CC=O>>CC(O)CC=O", "reaction_id":"r000004"},
                {"smiles":"CC=O.CC(O)CC=O>>CC(O)C(C=O)C(C)O", "reaction_id":"r000006"}

            ], seeds=["CC=O"],
            hits=["CC(O)C(C=O)C(O)C"]
        )
        self.assertEqual(results['r000004']['active'], True)
        self.assertEqual(results['r000006']['active'], True)
        
    def test_lp_prune_4(self):

        status, results = lp_prune(
            reactions=[
                {"smiles": "CC=O.CC=O>>CC(O)CC=O", "reaction_id":"r00"},
                {"smiles": "C=O.C=O>>CCCO", "reaction_id":"r01"}
            ], seeds=["C=O", "C=O"],
            hits=["CCCO"]
        )
        self.assertEqual(status, 1)
        self.assertEqual(results['r01']['active'], True)

    def test_lp_prune_5(self):

        status, results = lp_prune(
            reactions=[
                {"smiles": "CC=O.CC=O>>CC(O)CC=O", "reaction_id":"r000000"},
                {"smiles":"C=O.C=O>>CCO", "reaction_id":"r000001"}
            ], seeds=["N", "S"],
            hits=["CCO"]
        )
        self.assertEqual(status, -1)
    

    def test_prune_real_network_1(self):

        with open('tests/expansion-1.json') as f:
            u = json.load(f)

        status, results = lp_prune(
            reactions=u['reactions'],
            seeds=[
                "CC(=O)O", "C(=O)O", "CO", "N", "NO", "OS(=O)(=O)O", "OP(=O)(O)O", "C(=O)(O)O", "S"
            ], hits=["O=COC=O", "OC(O)S"]
        )
        self.assertEqual(status, 1)

    def test_prune_real_network_2(self):

        with open('tests/expansion-all.json') as f:
            u = json.load(f)

        status, results = lp_prune(
            reactions=u['reactions'],
            seeds=[
                "CC(=O)O", "C(=O)O", "CO", "N", "NO", "OS(=O)(=O)O", "OP(=O)(O)O", "C(=O)(O)O", "S"
            ], hits=["CC(O)(CC(=O)O)O=C(O)C(C(O)(O)O)C(O)(O)O"]
        )
        self.assertEqual(status, 1)

    def test_filter_hits(self):

        hits = filter_hits(
            reactions=[{"smiles":"C(O)=O>>CC"}, {"smiles":"CC>>CCC"}, {"smiles":"CCC>>CCN"}, {"smiles":"CCC>>CCCC"}],
            hits=['N', 'CCC']
        )
        self.assertEqual(len(hits), 1)

    def test_grow_and_prune(self):

        def test_function(self, rounds):
            reaction_rules = open('experiments/exp-2024-08-23/s001/rules.csv')
            seeds = open('experiments/exp-2024-08-23/s001/seeds.csv')
            
            reaction_rules = reaction_rules.read().split('\n')[:-1]
            seeds = seeds.read().split('\n')[:-1]
            reaction_rules = [dict(name='r{:06d}'.format(i), smarts=r) for i, r in enumerate(reaction_rules)]
            expansion_results = iterative_probabilistic_expansion(seeds=seeds, reaction_rules=reaction_rules, iterations=100, rounds=rounds)
            
            hits = expansion_results['seeds'][-20:]
            pruned_hits = filter_hits(reactions=expansion_results['reactions'], hits=hits)
            pruned_hits = [p['smiles'] for p in pruned_hits]
            status, sol = lp_prune(reactions=expansion_results['reactions'], seeds=seeds, hits=pruned_hits)
            self.assertEqual(status, 1)

        # test_function(self, 1)
        test_function(self, 2)
        test_function(self, 3)
        test_function(self, 4)



if __name__ == "__main__":
    unittest.main()