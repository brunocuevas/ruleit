import unittest
from ruleit.reaction_factory import _expansion


class TestFactory(unittest.TestCase):

    def test_factory(self):
        f = open('tests/seeds.txt')
        g = open('reaction-rules.yaml')
        _expansion(f, g, 'out')


if __name__ == '__main__':
    unittest.main()