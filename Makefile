all: test

test:
	python -m unittest discover -v

dep:
	python -m pip install -r requirements.txt

install: dep
	python -m pip install -e .

export:
	python -c 'from welib.tools.repo import *; export_figs_rec("welib/");'
