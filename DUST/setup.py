from setuptools import setup, Extension

module = Extension('dust_module',
                   sources=['dust_module.c', 'sdust.c'],
                   libraries=['z'])  

setup(name='dust_module',
      version='1.0',
      description='DUST module',
      ext_modules=[module])