VERSION = 1.0.0
REPO = variant-filtration-tools
MODULE = gdc_filtration_tools
BRANCH_NAME?=unknown

GIT_SHORT_HASH:=$(shell git rev-parse --short HEAD)

LONG_VERSION:=$(shell python3 setup.py -q capture_version --semver ${VERSION} --branch ${BRANCH_NAME})
PYPI_VERSION:=$(shell python3 setup.py -q print_version --pypi)
COMMIT_HASH:=$(shell python3 setup.py -q print_version --hash)

DOCKER_REPO := quay.io/ncigdc
DOCKER_IMAGE := ${DOCKER_REPO}/${REPO}:${LONG_VERSION}
DOCKER_IMAGE_COMMIT := ${DOCKER_REPO}/${REPO}:${COMMIT_HASH}
DOCKER_IMAGE_LATEST := ${DOCKER_REPO}/${REPO}:latest

.PHONY: version version-*
version:
	@echo --- VERSION: ${LONG_VERSION} ---

version-short:
	@echo ${VERSION}

version-long:
	@echo ${LONG_VERSION}

version-pypi:
	@echo ${PYPI_VERSION}

version-docker:
	@echo ${DOCKER_IMAGE}
	@echo ${DOCKER_IMAGE_COMMIT}
	@echo ${DOCKER_IMAGE_LATEST}

.PHONY: build build-* clean init init-* lint requirements run version
init: init-pip init-hooks

init-pip:
	@echo
	@echo -- Installing pip packages --
	pip3 install --no-cache-dir -r requirements.txt

init-hooks:
	@echo
	@echo -- Installing Precommit Hooks --
	pre-commit install

init-venv:
	@echo
	PIP_REQUIRE_VIRTUALENV=true pip3 install --upgrade pip-tools

clean:
	rm -rf ./build/
	rm -rf ./dist/
	rm -rf ./*.egg-info/

lint:
	@echo
	@echo -- Lint --
	python3 -m flake8 \
		--ignore=E501,F401,E302,E502,E126,E731,W503,W605,F841,C901 \
		${MODULE}/

run:
	bin/run

requirements: init-venv
	python3 setup.py -q capture_requirements --dev
	pip-compile -o requirements.txt requirements.in

.PHONY: build build-*

build: build-docker

build-docker:
	@echo
	@echo -- Building docker --
	python3 setup.py build
	mkdir -p dist
	cp -r build/lib/* dist/
	cp -r bin/ dist/
	cp -r tests/ dist/
	cp -f Makefile requirements.txt README.md setup.py dist/
	docker build . \
		--file ./Dockerfile \
		-t "${DOCKER_IMAGE_COMMIT}" \
		-t "${DOCKER_IMAGE}" \
		-t "${DOCKER_IMAGE_LATEST}"

.PHONY: test test-*
test: lint test-unit

test-unit:
	@echo
	@echo -- Unit Test --
	python3 -m pytest --cov-report term-missing --cov=${MODULE} tests/

test-docker:
	@echo
	@echo -- Running Docker Test --
	docker run --rm ${DOCKER_IMAGE_LATEST} test
