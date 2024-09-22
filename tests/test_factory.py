import unittest
from ruleit.expansion import _expansion, probablistic_expansion, safety_check, iterative_probabilistic_expansion
from rdkit.Chem import rdChemReactions as rdr
from yaml import Loader, Dumper, load, dump
import numpy as np


class TestFactory(unittest.TestCase):

    def test_factory(self):
        f = list(map(lambda x: x.strip(), open('tests/seeds.txt').readlines()))
        g = load(open('tests/reaction-rules.yaml'), Loader)
        u = _expansion(f, g, 100)
        len(u['discovered-reactions'])

    def test_factory_multiple_steps(self):
        f = list(map(lambda x: x.strip(), open('tests/seeds.txt').readlines()))
        g = load(open('tests/reaction-rules.yaml'), Loader)
        for i in range(4):
            u = _expansion(f, g, 1000)
            f = u['discovered-molecules']
        self.assertEqual(len(u['discovered-reactions']), 1000)

    def test_all_rules(self):
        f = list(map(lambda x: x.strip(), open('tests/seeds.test.txt').readlines()))
        
        g = list(map(lambda x: x.strip(), open('tests/reactions.test.txt').readlines()))
        g = [dict(smarts=item, name='r{:06d}'.format(i)) for i, item in enumerate(g)]
        g = dict(reactions=g)
        for i in range(3):
            
            u = _expansion(f, g, 10000)
            f = u['discovered-molecules']
        # self.assertEqual(len(u['discovered-reactions']), 1000)

    def test_probabilistic_expansion(self):
        f = list(map(lambda x: x.strip(), open('tests/seeds.txt').readlines()))
        g = load(open('tests/reaction-rules.yaml'), Loader)
        p = np.array([0.25, 0.25, 0.25, 0.25])
        
        u = probablistic_expansion(f, g, rule_probability=p, iterations=1000)
        px = len(list(filter(lambda x: x['rule'] == 'Isomerization', u['discovered-reactions']))) / len(u['discovered-reactions'])
        print("---")

        u = probablistic_expansion(f, g, rule_probability=p, iterations=10000)
        px = len(list(filter(lambda x: x['rule'] == 'Isomerization', u['discovered-reactions']))) / len(u['discovered-reactions'])
        print("---")

    def test_check_sanity(self):
        test_set = [
            (rdr.ReactionFromSmarts('CCC.N>>CCCC(N)O'), True),
            (rdr.ReactionFromSmarts('CCC.N>>CCCCCCC(C)(=N)O'), False),
            
            (rdr.ReactionFromSmarts('CCC.N>>CCCCCC(C)(C)C'), True),
            (rdr.ReactionFromSmarts('CCC.N>>BrCC(Br)CCCCCC(N)CCCC(Cl)CCCCCCCCCCCC'), False),
            (rdr.ReactionFromSmarts('CCC.N>>N#CCCCCCC(C)O'), True),
        ]
        conditions = dict(
            mass = 500.0,
            valence = {
                "C": [4]
            }
        )
        for r in test_set:
            self.assertEqual(safety_check(r[0], conditions), r[1])

    def test_probabilistic_capped(self):
        f = list(map(lambda x: x.strip(), open('tests/seeds.test.txt').readlines()))
        
        g = list(map(lambda x: x.strip(), open('tests/reactions.test.txt').readlines()))
        g = [dict(name='r{:06d}'.format(i), smarts=r) for i, r in enumerate(g)]
        #g = dict(reactions=g)
        conditions = dict(
            mass=500.0,
            valence=dict(
                C=[4],
                N=[3],
                O=[2]
            )
        )

        out = iterative_probabilistic_expansion(
            seeds=f, reaction_rules=g, iterations=1000, rounds=4,
            conditions=conditions
        )
        self.assertGreater(len(out['reactions']), 1)
        print(len(out['reactions']))

    def test_example_1(self):
        f = list(map(lambda x: x.strip(), open('data/example.1.seeds.csv').readlines()))
        g = list(map(lambda x: x.strip(), open('data/example.1.reactions.csv').readlines()))
        g = [dict(name='r{:06d}'.format(i), smarts=r) for i, r in enumerate(g)]
        #g = dict(reactions=g)
        conditions = dict(
            mass=500.0,
            valence=dict(
                C=[4],
                N=[3],
                O=[2]
            )
        )

        out = iterative_probabilistic_expansion(
            seeds=f, reaction_rules=g, iterations=10000, rounds=4,
            conditions=conditions
        )
        self.assertGreater(len(out['reactions']), 1)
        print(len(out['reactions']))



if __name__ == '__main__':
    unittest.main()