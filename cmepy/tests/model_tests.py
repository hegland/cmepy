import unittest
from test import test_support


# models will automatically validate once instantiated,
# so it should be sufficient to instantiate them for testing
# purposes

class ModelTests(unittest.TestCase):
    def test_burr08_models(self):
        from cmepy.models import burr08
        m = burr08.create_model_competing_clonotypes()

    def test_dsmts_models(self):
        from cmepy.models import dsmts
        m = dsmts.DSMTS_001_01
    
    def test_dual_enzymatic_models(self):
        from cmepy.models import dual_enzymatic
        m = dual_enzymatic.create_model()
        states = set(x for x in dual_enzymatic.gen_states())
    
    def test_gou07_models(self):
        from cmepy.models import gou07
        m1 = gou07.create_model_uni_dim()
        m2 = gou07.create_model_quad_autocat()
        m3 = gou07.create_model_quad_autocat(fixed_s = False)
    
    def test_michaelis_menten_models(self):
        from cmepy.models import michaelis_menten
        m = michaelis_menten.create_model_michaelis_menten()
    
    def test_mono_molecular_models(self):
        from cmepy.models import mono_molecular
        m1 = mono_molecular.A2B2A
        m2 = mono_molecular.A2B2C
    
    def test_munk08_models(self):
        from cmepy.models import munk08
        m = munk08.create_model_gene_toggle()

def suite():
    suite = unittest.TestLoader().loadTestsFromTestCase(ModelTests)
    return suite

def main():
    test_support.run_unittest(ModelTests)

if __name__ == '__main__':
    main()