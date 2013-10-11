from distutils.core import setup
#    url='http://pypi.python.org/pypi/WhiteMatterAnalysis/',
#    scripts=['bin/test1.py','bin/test2.py'],

setup(
    name='WhiteMatterAnalysis',
    version='0.1.0',
    author='Lauren O\'Donnell',
    author_email='odonnell@bwh.harvard.edu',
    packages=['whitematteranalysis', 'whitematteranalysis.test'],
    license='LICENSE.txt',
    description='Processing of whole-brain streamline tractography.',
    long_description=open('README.txt').read(),
)

