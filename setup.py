import os, glob
from setuptools import setup


# Functions read() was copied from Pip package.
def read(rel_path):
    here = os.path.abspath(os.path.dirname(__file__))
    # intentionally *not* adding an encoding option to open, See:
    #   https://github.com/pypa/virtualenv/issues/201#issuecomment-3145690
    with open(os.path.join(here, rel_path), 'r') as fp:
        return fp.read()


setup(
    name='prefactor',
    version='3.0',
    description='Prefactor: The LOFAR pre-facet calibration pipeline.',
    long_description=read("README.md"),
    long_description_content_type='text/markdown',
    url='https://github.com/lofar-astron/prefactor',
    license='GNU GPL 3',
    classifiers=[
        'Programming Language :: Python :: 3',
        'Development Status :: 3 - Alpha',
        'Natural Language :: English',
        'Intended Audience :: Science/Research',
        'Operating System :: POSIX :: Linux',
        'Topic :: Scientific/Engineering :: Astronomy'],
    platforms='any',
    install_requires=[
        'aplpy', 'astropy', 'bdsf', 'losoto', 'lsmtool',
        'matplotlib', 'numpy', 'python-casacore', 'scipy'],
    scripts=glob.glob('scripts/*'),
    data_files=[
        ('share/prefactor/rfistrategies', glob.glob('rfistrategies/*')),
        ('share/prefactor/skymodels', glob.glob('skymodels/*')),
        ('share/prefactor/solutions', glob.glob('solutions/*'))]
)
