FROM quay.io/jeremiahsavage/cdis_base

USER root

RUN apt-get update && apt-get install -y --force-yes \
    wget \
    python-dev \
    unzip \
    cmake \
    libncurses-dev \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev

RUN apt-get clean \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

USER ubuntu

ENV HOME /home/ubuntu

RUN mkdir -p ${HOME}/tools/

ENV variant-filtration-tool 2.0 

WORKDIR ${HOME}/tools/

## Install HTSLIB
ENV VERSION 1.3.2
ENV NAME htslib
ENV URL "https://github.com/samtools/htslib/releases/download/1.3.2/htslib-1.3.2.tar.bz2"

RUN wget -q -O - ${URL} | tar -xjf - && \
    cd ${NAME}-${VERSION} && \
    ./configure && \
    make

ENV PATH ${PATH}:${HOME}/tools/${NAME}-${VERSION}/

## Install bam-readcount
RUN wget https://github.com/genome/bam-readcount/archive/v0.7.4.tar.gz && \
    tar -xzf v0.7.4.tar.gz && \
    rm -f v0.7.4.tar.gz && \
    mv bam-readcount-0.7.4 ${HOME}/tools/bam-readcount

RUN cd ${HOME}/tools/bam-readcount && mkdir build

USER root
RUN cd ${HOME}/tools/bam-readcount/build && \
    cmake ../ && \
    make deps && \
    make -j && \
    make install

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
