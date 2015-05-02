from distutils.core import setup
#    url='http://pypi.python.org/pypi/WhiteMatterAnalysis/',
#    scripts=['bin/test1.py','bin/test2.py'],
from Cython.Build import cythonize
import numpy

setup(
    name='WhiteMatterAnalysis',
    version='0.1.0',
    author='Lauren O\'Donnell',
    author_email='odonnell@bwh.harvard.edu',
    packages=['whitematteranalysis'],
    license='LICENSE.txt',
    description='Processing of whole-brain streamline tractography.',
    long_description=open('README.md').read(),
    ext_modules = cythonize("whitematteranalysis/*.pyx"),
    include_dirs=[numpy.get_include()]
)

