all: test

test:
	python -m unittest discover


install:
	python -m pip install -r requirements.txt
