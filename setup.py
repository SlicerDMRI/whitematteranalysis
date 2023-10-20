import glob
import os
from itertools import chain
from setuptools import find_packages, setup

#    url='http://pypi.python.org/pypi/whitematteranalysis/',
#    scripts=['bin/test1.py','bin/test2.py'],

# ext_modules work-around: https://stackoverflow.com/a/38057196
# get_include work-around: https://stackoverflow.com/a/21621689

import sys
if sys.platform == 'win32':
    # force setuptools to use the MSVC compiler environment
    # otherwise it will check version for VC9, and error out
    os.environ['MSSdk'] = '1'
    os.environ['DISTUTILS_USE_SDK'] = '1'


with open("requirements.txt") as f:
    required_dependencies = f.read().splitlines()
    external_dependencies = []
    for dependency in required_dependencies:
        if dependency[0:2] == "-e":
            repo_name = dependency.split("=")[-1]
            repo_url = dependency[3:]
            external_dependencies.append("f{repo_name} @ {repo_url}")
        else:
            external_dependencies.append(dependency)


setup(
    name='whitematteranalysis',
    version='0.3.0',
    author='Fan Zhang and Lauren O\'Donnell',
    author_email='odonnell@bwh.harvard.edu',
    packages=find_packages(),
    license='LICENSE.txt',
    description='Processing of whole-brain streamline tractography.',
    long_description=open('README.md').read(),
    long_description_content_type="text/markdown",
    install_requires=external_dependencies,
    extras_require={'doc': ['sphinx', 'sphinx-argparse', 'sphinx_rtd_theme']},
    scripts=list(chain.from_iterable([
        glob.glob("bin/[a-zA-Z]*.py"),
        glob.glob("utilities/[a-zA-Z]*.py"),
        ]))
)
