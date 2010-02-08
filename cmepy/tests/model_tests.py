import unittest
from test import test_support

from cmepy import validate

class ModelTests(unittest.TestCase):
    def setUp(self):
        from cmepy.models.gou07 import GOU07_A, GOU07_B, GOU07_C
        from cmepy.models.munk08 import MUNK08_A, MUNK08
        from cmepy.models.michaelis_menten import MM_SIMPLE
        from cmepy.models.dsmts import DSMTS_001_01
        from cmepy.models.mono_molecular import A2B2C, A2B2A
        from cmepy.models.burr08_model import BURR08
        self.models = [
            GOU07_A, GOU07_B, GOU07_C,
            MUNK08_A, MUNK08,
            MM_SIMPLE,
            DSMTS_001_01,
            A2B2C, A2B2A,
            BURR08,
        ]
    
    def test_model_validity(self):
        for m in self.models:
            validate.model(m)

def suite():
    suite = unittest.TestLoader().loadTestsFromTestCase(ModelTests)
    return suite

def main():
    test_support.run_unittest(ModelTests)

if __name__ == '__main__':
    main()
