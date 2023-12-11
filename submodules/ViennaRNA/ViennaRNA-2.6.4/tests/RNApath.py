import os, sys

VERSION_NUMBER = "2.6.4"

def addSwigInterfacePath(version=3):
    if version == 2:
        sys.path.insert(0, os.path.join('..', 'interfaces/Python2'))
    else:
        sys.path.insert(0, os.path.join('..', 'interfaces/Python'))

def getDataDirPath():
    return os.path.join("..", "tests/data/")
