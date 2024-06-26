import unittest
from ruleit.expansion import _expansion
from yaml import Loader, Dumper, load, dump


class TestFactory(unittest.TestCase):

    def test_factory(self):
        f = list(map(lambda x: x.strip(), open('tests/seeds.txt').readlines()))
        g = load(open('tests/reaction-rules.yaml'), Loader)
        u = _expansion(f, g, 'out')
        len(u['discovered-reactions'])

    def test_factory_multiple_steps(self):
        f = list(map(lambda x: x.strip(), open('tests/seeds.txt').readlines()))
        g = load(open('tests/reaction-rules.yaml'), Loader)
        for i in range(4):
            u = _expansion(f, g, 1000)
            f = u['discovered-molecules']
        self.assertEqual(len(u['discovered-reactions']), 1000)


if __name__ == '__main__':
    unittest.main()