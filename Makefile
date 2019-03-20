all: test

test:
	python -m unittest discover

dep:
	python -m pip install -r requirements.txt

install: dep
	python -m pip install -e .
