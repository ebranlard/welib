all: test

test:
	python -m unittest discover -v

dep:
	python -m pip install -r requirements.txt

install: dep
	python -m pip install -e .

export:
	python -c 'from welib.tools.repo import *; export_figs_rec("welib/");'



# --- Tools for pypi
clean:
	rm -rf __pycache__
	rm -rf *.egg-info
	rm -rf build
	rm -rf dist

pypi0:
	python setup.py sdist bdist_wheel

pypi1:
	twine upload dist/*

pipun:
	pip uninstall -y -r requirements.txt

example:
	python welib/beams/examples/Ex1_BeamModes.py
