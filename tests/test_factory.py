import unittest
from ruleit.expansion import _expansion
from yaml import Loader, Dumper, load, dump


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



if __name__ == '__main__':
    unittest.main()