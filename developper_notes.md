# Developer Notes 

## Release Steps
### --- Main Steps (create a package and uploat to pypi): 
change VERSION in setup.py 
```bash
python setup.py sdist
twine upload dist/*
```

### --- Detailed steps
Go to folder
```bash
cd path/to/welib
```

Create a source distribution
```bash
python setup.py sdist
Install twine
```bash
pip install twine
```
Run twine to upload to Pypi (will ask for username and password)
```bash
twine upload dist/*
```



## After clone / first time
Add `.gitconfig` to your path, to apply filters on jupyter notebooks
```bash
git config --local include.path ../.gitconfig
```



## Conda 

conda-forge, 
make a pull request there.
 - setup build script (https://conda-forge.org/docs/maintainer/adding_pkgs.html). 
    see e.g. FLORIS (https://github.com/conda-forge/floris-feedstock/blob/master/recipe/meta.yaml).
 - make a pull request to https://github.com/conda-forge/staged-recipes/ 
