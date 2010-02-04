import unittest
from test import test_support


# models will automatically validate once instantiated,
# so it should be sufficient to instantiate them for testing
# purposes

class ModelTests(unittest.TestCase):
    def test_burr08_models(self):
        from cmepy.models import burr08
        m = burr08.create_model_competing_clonotypes()
    
    def test_mono_molecular_models(self):
        from cmepy.models import mono_molecular
        m1 = mono_molecular.A2B2A
        m2 = mono_molecular.A2B2C
    
    def test_gou07_models(self):
        from cmepy.models import gou07
        m1 = gou07.create_model_uni_dim()
        m2 = gou07.create_model_quad_autocat()
        m3 = gou07.create_model_quad_autocat(fixed_s = False)

def suite():
    suite = unittest.TestLoader().loadTestsFromTestCase(ModelTests)
    return suite

def main():
    test_support.run_unittest(ModelTests)

if __name__ == '__main__':
    main()