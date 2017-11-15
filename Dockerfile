FROM ubuntu:16.04

MAINTAINER Kyle Hernandez <kmhernan@uchicago.edu>

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get -y update && apt-get install -y --force-yes \
    build-essential \
    autoconf \
    zlibc zlib1g-dev \
    libncurses-dev \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    python3.5 \ 
    python3.5-dev \
    python3-pip \
    wget \
    time \
    git-core

RUN apt-get clean \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

ENV variant-filtration-tool 2.5

## Install vt
RUN wget -q -O - https://github.com/atks/vt/archive/0.5772.tar.gz  
## Install variant-filtration-tool
WORKDIR ${HOME}
RUN mkdir -p ${HOME}/tools/variant-filtration-tool
ADD utils ${HOME}/tools/variant-filtration-tool/
ADD requirements.txt ${HOME}/tools/variant-filtration-tool/

RUN /bin/bash -c "source ${HOME}/.local/bin/virtualenvwrapper.sh \
    && mkvirtualenv --python=/usr/bin/python2.7 p2 \
    && cd ~/tools/variant-filtration-tool \
    && pip install -r ./requirements.txt \
    && echo source ${HOME}/.virtualenvs/p2/bin/activate >> ${HOME}/.bashrc"

WORKDIR ${HOME}
