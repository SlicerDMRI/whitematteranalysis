#!/bin/sh

INSTALLFOLDER=/Applications/WhiteMatterAnalysis
WMAFOLDER=$(pwd)/..
PYTHONVERSION=3.6.0

mkdir -p $INSTALLFOLDER
cd $INSTALLFOLDER

git clone https://github.com/pyenv/pyenv
PREFIX=$INSTALLFOLDER/python-build pyenv/plugins/python-build/install.sh

$INSTALLFOLDER/python-build/bin/python-build $PYTHONVERSION $INSTALLFOLDER/python
$INSTALLFOLDER/python/bin/pip install $WMAFOLDER 


