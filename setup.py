from setuptools import setup, find_packages

VERSION='0.0.2'


EXTRAS = {
    'docs': {
        'readthedocs-sphinx-ext>=0.5.15',
        'Sphinx>=2.0',
        'sphinxcontrib-napoleon>=0.7'
    }
}

setup(
    name='welib',
    version=VERSION,
    description='Wind Energy Library',
    long_description="""
Wind energy library: suite of python tools for aero-servo-hydro-elasticity (aerodynanmics, controls, hydrodynamics, structure/elasticity) and wind energy.""",
    long_description_content_type = 'text/markdown',
    author='Emmanuel Branlard',
    author_email='lastname@gmail.com',
    url='http://github.com/ebranlard/welib/',
    license='MIT',
    python_requires=">=3.6",
    packages=find_packages(exclude=["tests"]),
    install_requires=[
        'matplotlib', 
        'xlrd',
        'numpy',
        'pandas', 
        'future', 
        'chardet',
        'scipy', 
        'sympy'
    ],
    extras_require       = EXTRAS,
    include_package_date = True,
    zip_safe=False,
    classifiers=[
              'Development Status :: 5 - Production/Stable',
              'Environment :: Console',
              'Intended Audience :: Science/Research',
              'Intended Audience :: Education',
              'Intended Audience :: End Users/Desktop',
              'Intended Audience :: Developers',
              'License :: OSI Approved :: MIT License',
              'Operating System :: MacOS :: MacOS X',
              'Operating System :: Microsoft :: Windows',
              'Operating System :: POSIX',
              'Operating System :: Unix',
              'Programming Language :: Python',
              'Topic :: Scientific/Engineering',
              'Topic :: Scientific/Engineering :: Atmospheric Science',
              'Topic :: Scientific/Engineering :: Hydrology',
              'Topic :: Scientific/Engineering :: Mathematics',
              'Topic :: Scientific/Engineering :: Physics'
              ],
)
