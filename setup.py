#!python

"""
CmePy setup script
"""


# ensure setuptools is present
# n.b. no_fake flag should prevent this script
#      from 'patching' the currently installed
#      setuptools, if any

from distribute_setup import use_setuptools
use_setuptools(no_fake = True)


VERSION = '0.2.0'


from setuptools import setup, find_packages
setup(
    name = 'cmepy',
    version = VERSION,
    
    package_data = {'':['*.txt'], 'doc':['*.pdf']},
    
    test_suite = "test_all",

    packages = find_packages(),
)

