#### Make configuration:

## Use Bash as default shell, and in strict mode:
SHELL := /bin/bash
.SHELLFLAGS = -ec

## If the parent env doesn't ste TMPDIR, do it ourselves:
TMPDIR ?= /tmp

## Users can override this variable from the command line,
## to install MP binaries somewhere other than /usr/local,
## if they lack root privileges:
DESTDIR ?= /usr/local

## This makes all recipe lines execute within a shared shell process:
## https://www.gnu.org/software/make/manual/html_node/One-Shell.html#One-Shell
.ONESHELL:

## If a recipe contains an error, delete the target:
## https://www.gnu.org/software/make/manual/html_node/Special-Targets.html#Special-Targets
.DELETE_ON_ERROR:

## This is necessary to make sure that these intermediate files aren't clobbered:
.SECONDARY:


### Local Definitions:
PYTHON ?= python3

### Install Singularity

OS := linux
ARCH := amd64
GO-VERSION := 1.17.6
SY-VERSION := 3.9.3
singularity-install:
	sudo apt-get update
	sudo apt-get install -y \
	   build-essential \
	   libseccomp-dev \
	   pkg-config \
	   squashfs-tools \
	   cryptsetup
	cd $(TMPDIR)
	wget https://dl.google.com/go/go$(GO-VERSION).$(OS)-$(ARCH).tar.gz # Downloads the required Go package
	sudo tar -C /usr/local -xzvf go$(GO-VERSION).$(OS)-$(ARCH).tar.gz  # Extracts the archive
	rm go$(GO-VERSION).$(OS)-$(ARCH).tar.gz                            # Deletes the ``tar`` file
	export PATH=$$PATH:/usr/local/go/bin
	wget https://github.com/sylabs/singularity/releases/download/v${SY-VERSION}/singularity-ce-${SY-VERSION}.tar.gz
	tar -xzf singularity-ce-${SY-VERSION}.tar.gz
	rm singularity-ce-${SY-VERSION}.tar.gz
	cd singularity-ce-${SY-VERSION}
	./mconfig
	make -C builddir
	sudo make -C builddir install


### Container Automation
docker-start:
	sudo systemctl start docker

## If git_branch is empty, it probably means that we are building off of a tagged version.
## So, we then grab the tag string:
docker-build:
	git_branch=$$(git symbolic-ref --short -q HEAD) \
		|| git_branch=$$(git describe --tags)
	sudo docker build --network=host \
			--build-arg git_branch=$$git_branch \
			-t quay.io/hallamlab/metapathways:$$git_branch -f docker/Dockerfile .

docker-run:
	git_branch=$$(git symbolic-ref --short -q HEAD)
	sudo docker run -it --network=host --rm \
		-v $(CURDIR):/input \
		-v $(CURDIR)/out:/output \
		quay.io/hallamlab/metapathways:$$git_branch /bin/bash

docker-test:
	git_branch=$$(git symbolic-ref --short -q HEAD)
	sudo docker run --rm --workdir /tmp \
		quay.io/hallamlab/metapathways:$$git_branch \
		/bin/bash -c 'metapathways build_db --test && metapathways run --test'


docker-deploy:
	sudo docker login quay.io
	sudo docker push quay.io/hallamlab/metapathways:dev

docker-fetch:
	sudo docker pull quay.io/hallamlab/metapathways

singularity-local-build:
	git_branch=$$(git symbolic-ref --short -q HEAD) \
		|| git_branch=$$(git describe --tags)
	sudo singularity build metapathways-$$git_branch.sif docker-daemon://quay.io/hallamlab/metapathways:$$git_branch

singularity-local-shell:
	singularity shell metapathways-dev.sif

singularity-docker-build:
	sudo singularity build metapathways-dev.sif docker://quay.io/hallamlab/metapathways:dev

singularity-docker-shell:
	singularity shell docker://quay.io/hallamlab/metapathways:dev



### Packaging:
##

### PyPi:

create-package: clean-package
	$(PYTHON) -m pip install --user --upgrade setuptools wheel twine
	$(PYTHON) setup.py sdist bdist_wheel --universal
	cp -r dist docker/dist

clean-package:
	rm -rf dist metapathways.egg-info build docker/dist

install-package:
	$(PYTHON) -m pip install --user .

install-dev-package:
	$(PYTHON) -m pip install --user --upgrade -e .

install-dist-package:
	$(PYTHON) -m pip install --user dist/metapathways-*-py*-none-any.whl

deploy-package-to-pypi:
	twine upload dist/*

### Conda:

CONDA ?= mamba
ENV_NAME ?= mpw_dev
CONDA_RUN = $(CONDA) run -n $(ENV_NAME)

conda-deploy:
	$(PYTHON) scripts/release.py upload-conda

clean-dev-env:
	@echo ">>> Removing conda env '$(ENV_NAME)' (if it exists)"
	CONDA_NO_PLUGINS=true $(CONDA) env remove -n $(ENV_NAME) -y || true

create-dev-env:
	@echo ">>> Ensuring conda env '$(ENV_NAME)' exists (using $(CONDA))"
	# Increase open-file limit to avoid 'Too many open files' (Errno 24)
	@ulimit -n 65535 || true
	@if $(CONDA) env list | awk '{print $$1}' | grep -qx "$(ENV_NAME)"; then \
		echo ">>> Environment '$(ENV_NAME)' already exists – updating from YAMLs"; \
		CONDA_NO_PLUGINS=true $(CONDA) env update -n $(ENV_NAME) --file ./docker/conda_base.yml; \
		CONDA_NO_PLUGINS=true $(CONDA) env update -n $(ENV_NAME) --file ./docker/conda_dev.yml; \
	else \
		echo ">>> Environment '$(ENV_NAME)' does not exist – creating"; \
		CONDA_NO_PLUGINS=true $(CONDA) env create --no-default-packages -n $(ENV_NAME) --file ./docker/conda_base.yml; \
		CONDA_NO_PLUGINS=true $(CONDA) env update -n $(ENV_NAME) --file ./docker/conda_dev.yml; \
	fi

# Build PyPI package using the env
conda-pip-package:
	$(CONDA_RUN) python -m pip install --upgrade setuptools wheel twine
	$(CONDA_RUN) python setup.py sdist bdist_wheel --universal
	cp -r dist docker/dist

# Build and integration-test the exact release package.
conda-build-package:
	$(PYTHON) scripts/release.py build

# Release controller: validates the exact package before upload.
.PHONY: release-prepare release-build release-publish release-upload-conda full-build
release-prepare:
	@test -n "$(VERSION)" || { echo "Use make release-prepare VERSION=3.5.0"; exit 1; }
	$(PYTHON) scripts/release.py prepare "$(VERSION)"

release-build:
	$(PYTHON) scripts/release.py build

release-publish:
	$(PYTHON) scripts/release.py publish

release-upload-conda:
	$(PYTHON) scripts/release.py upload-conda

# Kept as an alias; sequencing prevents uploads racing builds under make -j.
full-build: release-build
	$(PYTHON) scripts/release.py upload-conda

# Legacy aliases use the same validation gates as the release controller.
create-conda: release-build

deploy-conda: release-upload-conda

### Docs:

docs-local:
	sphinx-build ./docs/src ./docs/build

### Build & Install Extensions
##
##

extensions-build:
	$(MAKE) -C extensions clean
	$(MAKE) -C extensions

## Deprecated, see Python package.
extensions-build-install: extensions-build
	mkdir -p $(DESTDIR)/bin
	cp extensions/FAST/fast*            $(DESTDIR)/bin
	cp extensions/metacount/metacount   $(DESTDIR)/bin

## Deprecated, these binaries need to get revisited.
extensions-install:
	mkdir -p $(DESTDIR)/bin
	$(MAKE) -C extensions/metacount clean
	$(MAKE) -C extensions/metacount 
	cp extensions/FAST/fastal $(DESTDIR)/bin/fastal
	cp extensions/FAST/fastdb $(DESTDIR)/bin/fastdb
	cp extensions/metacount/metacount $(DESTDIR)/bin/metacount
	chmod 755 $(DESTDIR)/bin/fastal
	chmod 755 $(DESTDIR)/bin/fastdb
	chmod 755 $(DESTDIR)/bin/metacount

# Validated release containers (Docker plus Apptainer).
.PHONY: release-containers release-quay
release-containers:
	$(PYTHON) scripts/release.py container-build

release-quay:
	$(PYTHON) scripts/release.py container-push
