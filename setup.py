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
    python_requires='>=3.5',
    install_requires=['numpy', 'scipy', 'PyQt5', 'pyqtgraph', 'lmfit',
                      "importlib_resources ; python_version<'3.7'"],
    package_data={'LiquidDiffract': ['resources/*', 'resources/icons/*',
                                     'resources/hubbel_compton/*',
                                     'resources/docs/*',
                                     'scripts/*']},
    zip_safe=True,
    entry_points={'console_scripts':
                  ['LiquidDiffract=LiquidDiffract.LiquidDiffract:main']}
)
