import unittest






class Testgff(unittest.TestCase):

    """Docstring for Testgff. """

    def SetUp(self):
        """TODO: to be defined1. """
        unittest.Test.__init__(self)

        
    def test_parser(self):
        parser = parse_args(['-g', '-ra', '-s', '-l','-o'])
        self.assertTrue(parser.long)

