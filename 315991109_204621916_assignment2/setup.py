from setuptools import setup, Extension

setup(name='mykmeanssp',
    version='0.0',
    description='kmeans exercise', 
    ext_modules=[Extension('mykmeanssp', sources=['kmeans.c'])])
