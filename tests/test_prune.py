from ruleit.lp_prune import *
import unittest


class TestLP(unittest.TestCase):

    def test_lp_prune_1(self):

        status, results = formalize_problem(
            reactions=[
                {"smiles":"C>>CC"},
                {"smiles":"CC>>CCC"},
                {"smiles":"CCC>>CCN"},
                {"smiles":"CCC>>CCCC"}
            ], seeds=["C"],
            hits=["CCN"]
        )
        self.assertEqual(results['r000000']['active'], True)
        self.assertEqual(results['r000001']['active'], True)
        self.assertEqual(results['r000002']['active'], True)
        self.assertEqual(results['r000003']['active'], False)

    def test_lp_prune_2(self):

        status, results = formalize_problem(
            reactions=[
                {"smiles":"C(O)=O>>CC"},
                {"smiles":"CC>>CCC"},
                {"smiles":"CCC>>CCN"},
                {"smiles":"CCC>>CCCC"}

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

        status, results = formalize_problem(
            reactions=[
                {"smiles": "CC=O.CC=O>>CC(O)CC=O"},
                {"smiles":"CC=O.CC(O)CC=O>>CC(O)C(C=O)C(C)O"}

            ], seeds=["CC=O"],
            hits=["CC(O)C(C=O)C(O)C"]
        )
        self.assertEqual(results['r000000']['active'], True)
        self.assertEqual(results['r000001']['active'], True)
        
    def test_lp_prune_4(self):

        status, results = formalize_problem(
            reactions=[
                {"smiles": "CC=O.CC=O>>CC(O)CC=O"},
                {"smiles": "C=O.C=O>>CCCO"}
            ], seeds=["C=O", "C=O"],
            hits=["CCCO"]
        )
        self.assertEqual(status, 1)
        self.assertEqual(results['r000001']['active'], True)

    def test_lp_prune_5(self):

        status, results = formalize_problem(
            reactions=[
                {"smiles": "CC=O.CC=O>>CC(O)CC=O"},
                {"smiles":"C=O.C=O>>CCO"}
            ], seeds=["N", "S"],
            hits=["CCO"]
        )
        self.assertEqual(status, -1)
        

if __name__ == "__main__":
    unittest.main()