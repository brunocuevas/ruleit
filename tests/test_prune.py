from ruleit.lp_prune import *
import unittest


class TestLP(unittest.TestCase):

    def test_lp_prune_1(self):

        status, results = formalize_problem(
            reactions=[
                {"smiles":"A>>B"},
                {"smiles":"B>>C"},
                {"smiles":"C>>D"},
                {"smiles":"C>>E"}
            ], seeds=["A"],
            hits=["D"]
        )
        self.assertEqual(results['r000000']['active'], True)
        self.assertEqual(results['r000001']['active'], True)
        self.assertEqual(results['r000002']['active'], True)
        self.assertEqual(results['r000003']['active'], False)

    def test_lp_prune_2(self):

        status, results = formalize_problem(
            reactions=[
                {"smiles":"A>>B"},
                {"smiles":"B>>C"},
                {"smiles":"C>>D"},
                {"smiles":"C>>E"},
                {"smiles":"E.A>>F"},
                {"smiles":"A.A>>G"},

            ], seeds=[],
            hits=["G"]
        )
        if status == 1:
            self.assertEqual(results['r000000']['active'], False)
            self.assertEqual(results['r000001']['active'], False)
            self.assertEqual(results['r000002']['active'], False)
            self.assertEqual(results['r000003']['active'], False)
            self.assertEqual(results['r000004']['active'], False)
            self.assertEqual(results['r000005']['active'], False)
        else:
            self.assertEqual(status, -1)

    def test_lp_prune_3(self):

        status, results = formalize_problem(
            reactions=[
                {"smiles":"A>>B"},
                {"smiles":"B>>C"},
                {"smiles":"C>>D"},
                {"smiles":"C>>E"},
                {"smiles":"E.A>>F"},
                {"smiles":"A.A>>G"},

            ], seeds=["A"],
            hits=["G", "F"]
        )
        self.assertEqual(results['r000000']['active'], True)
        self.assertEqual(results['r000001']['active'], True)
        self.assertEqual(results['r000002']['active'], False)
        self.assertEqual(results['r000003']['active'], True)
        self.assertEqual(results['r000004']['active'], True)
        self.assertEqual(results['r000005']['active'], True)

    def test_lp_prune_4(self):

        status, results = formalize_problem(
            reactions=[
                {"smiles": "CC=O.CC=O>>CC(O)CC=O"},
                {"smiles": "C=O.C=O>>C=O=CO"}
            ], seeds=["C=O", "C=O"],
            hits=["C=O=CO"]
        )
        self.assertEqual(status, 1)
        self.assertEqual(results['r000001']['active'], True)

    def test_lp_prune_5(self):

        status, results = formalize_problem(
            reactions=[
                {"smiles": "CC=O.CC=O>>CC(O)CC=O"},
                {"smiles":"C=O.C=O>>C=O=CO"}
            ], seeds=["N", "S"],
            hits=["C=O=CO"]
        )
        self.assertEqual(status, -1)
        

if __name__ == "__main__":
    unittest.main()