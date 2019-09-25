from setuptools import setup, find_packages
from os import path

here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='LiquidDiffract',
    version='0.1',
    description='GUI program to treat experimental X-ray diffraction data of liquid structures',
    long_description=long_description,
    long_description_content_type='text/markdown',
    license='GPL',
    author='Benedict J Heinen',
    author_email='benedict.heinen@gmail.com',
    url='https://github.com/bjheinen/LiquidDiffract',
    packages=find_packages(),
    python_requires='>=3.5',
    install_requires=['numpy','scipy','PyQt5'],
    package_data={'LiquidDiffract': ['data/*', 'data/icons/*', 'data/hubbel-compton/*']},
    entry_points={'console_scripts': ['LiquidDiffract=LiquidDiffract.LiquidDiffract:main']}
)



