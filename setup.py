from setuptools import setup, find_packages, Extension

k_module = Extension('mykmeanssp', sources=['kmeans.c'])
spk_module = Extension('spkModule', sources=['spkmeansmodule.c'])

setup(
    name='spkmeans',
    version='1.0',
    author='Panker_Gorfung',
    description='spkmeans project', 
    packages=find_packages(),
    ext_modules=[spk_module, k_module]
    )
