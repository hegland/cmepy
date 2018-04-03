"""
unit tests for cmepy.models sub-package
"""

import unittest


# models will automatically validated once instantiated,
# so it should be sufficient to instantiate them for testing
# purposes

class ModelTests(unittest.TestCase):
    def test_burr08_models(self):
        from cmepy.models import burr08
        m = burr08.create_model()
        phi = burr08.create_time_dependencies()
    
    def test_catalytic_reaction_models(self):
        from cmepy.models import catalytic_reaction
        initial_counts = { 'A' : 20, 'D' : 15 }
        m = catalytic_reaction.create_model(**initial_counts)
        states = set(catalytic_reaction.gen_states(**initial_counts))
        
    def test_dsmts_models(self):
        from cmepy.models import dsmts
        m = dsmts.DSMTS_001_01
    
    def test_dual_enzymatic_models(self):
        from cmepy.models import dual_enzymatic
        m = dual_enzymatic.create_model()
        states = set(dual_enzymatic.gen_states())
    
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
    
    def test_pap_pili_models(self):
        from cmepy.models import pap_pili
        m = pap_pili.create_model()
        states = set(pap_pili.gen_states())
    
    def test_transcription_regulation_models(self):
        from cmepy.models import transcription_regulation
        m = transcription_regulation.create_model()
        phi = transcription_regulation.create_time_dependencies()

def suite():
    suite = unittest.TestLoader().loadTestsFromTestCase(ModelTests)
    return suite

def main():
    unittest.run(ModelTests)

if __name__ == '__main__':
    main()
