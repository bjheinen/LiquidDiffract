from setuptools import setup, find_packages
from os import path

here = path.abspath(path.dirname(__file__))

# Get the version number from version.py
version = {}
with open(path.join(here, 'LiquidDiffract/version.py')) as fp:
    exec(fp.read(), version)

# Get the long description from the README file
with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='LiquidDiffract',
    version=version['__version__'],
    description='Liquid structure/total X-ray scattering analysis',
    long_description=long_description,
    long_description_content_type='text/markdown',
    license='GPL',
    author='Benedict J Heinen',
    author_email='benedict.heinen@gmail.com',
    url='https://github.com/bjheinen/LiquidDiffract',
    packages=find_packages(),
    python_requires='>=3.8',
    install_requires=['numpy>=1.20.0', 'scipy>=1.6.0', 'PyQt6>=6.2.0',
                      'qtpy>=2.3.0', 'pyqtgraph>=0.13.0', 'lmfit>=1.0.2',
                      'packaging>=14.1', 'importlib_resources;python_version<"3.9"'],
    package_data={'LiquidDiffract': ['resources/*', 'resources/icons/*',
                                     'resources/hubbel_compton/*',
                                     'resources/docs/*',
                                     'scripts/*']},
    zip_safe=True,
    entry_points={'console_scripts':
                  ['LiquidDiffract=LiquidDiffract.LiquidDiffract:main']}
)
