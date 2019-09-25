#pip install git+https://github.com/bjheinen/LiquidDiffract









from setuptools import setup, find_namespace_packages

with open("README.md", 'r') as f:
    long_description = f.read()

setup(
    name='LiquidDiffract',
    version='0.1',
    description='GUI program to treat experimental X-ray diffraction data of liquid structures',
    long_description=long_description,
    license='GPL',
    author='Benedict J Heinen',
    author_email='benedict.heinen@gmail.com',
    url='https://github.com/bjheinen/LiquidDiffract',
    packages=find_namespace_packages(include=['LiquidDiffract.*']),
    install_requires=['numpy','scipy','PyQt5']

)
