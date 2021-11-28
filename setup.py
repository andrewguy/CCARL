import io
import os
import sys
from setuptools import setup
from setuptools.command.test import test as TestCommand

# Version is defined in ccarl/__init__.py
__version__ = "Undefined"
for line in open('ccarl/__init__.py'):
    if (line.startswith('__version__')):
        exec(line.strip())

here = os.path.abspath(os.path.dirname(__file__))


def read(*filenames, **kwargs):
    encoding = kwargs.get('encoding', 'utf-8')
    sep = kwargs.get('sep', '\n')
    buf = []
    for filename in filenames:
        with io.open(filename, encoding=encoding) as f:
            buf.append(f.read())
    return sep.join(buf)


LONG_DESCRIPTION = read('README.md', 'CHANGES.txt')


class PyTest(TestCommand):
    def finalize_options(self):
        TestCommand.finalize_options(self)
        self.test_args = []
        self.test_suite = True

    def run_tests(self):
        import pytest
        errcode = pytest.main(self.test_args)
        sys.exit(errcode)


setup(
    name='ccarl',
    version=__version__,
    url='https://github.com/andrewguy/CCARL',
    download_url='https://github.com/andrewguy/CCARL/archive/refs/tags/{0}.tar.gz'.format(__version__),
    author='Andrew Guy',
    tests_require=['pytest'],
    setup_requires=[],
    install_requires=['sklearn',
                      'pandas',
                      'numpy',
                      'networkx',
                      'matplotlib',
                      'statsmodels',
                      'pyxdameraulevenshtein',
                      'pyparsing',
                      'click'
                      ],
    cmdclass={'test': PyTest},
    author_email='andrewjguy42@gmail.com',
    description='A package for identification of glycan motifs',
    long_description=LONG_DESCRIPTION,
    packages=['ccarl'],
    include_package_data=True,
    platforms='any',
    classifiers=[
        'Programming Language :: Python',
        'Development Status :: 4 - Beta',
        'Natural Language :: English',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        ],
    extras_require={
        'testing': ['pytest'],
    },
    entry_points={
              'console_scripts': [
                  'ccarl-validate-cfg-structures = ccarl.cli:validate_structures',
                  'ccarl-identify-binders = ccarl.cli:identify_binders',
                  'ccarl-identify-motifs = ccarl.cli:identify_motifs',
                  'ccarl-predict-binding = ccarl.cli:predict_binding'
              ]
          }
)
