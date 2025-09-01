.PHONY: build test venv clean install upload

BIN=venv/bin/
PROJECT_NAME=cenplot

test:
	$(BIN)python3 -m pip install pytest
	$(BIN)python3 -m pytest -vv

build:
	$(BIN)python3 -m pip install --upgrade build
	$(BIN)python3 -m build

install:
	$(BIN)python3 -m pip uninstall -y $(PROJECT_NAME)
	$(BIN)python3 -m pip install $(shell find dist -name "*.whl" | sort | head -1) --no-input

venv:
	python3 -m venv venv

dev:
	$(MAKE) venv
	$(BIN)python3 -m pip install -r requirements.txt pytest pdoc pre-commit

clean:
	rm -rf dist/ venv/ .*cache/ *.egg-info/

upload:
	$(BIN)python3 -m pip install --upgrade twine
	$(BIN)python3 -m twine upload --repository pypi dist/*
