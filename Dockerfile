FROM quay.io/jeremiahsavage/cdis_base

USER root

RUN apt-get update && apt-get install -y --force-yes \
    wget \
    python-dev \
    unzip \
    cmake \
    libncurses-dev

RUN apt-get clean \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

USER ubuntu

ENV HOME /home/ubuntu

RUN mkdir -p ${HOME}/tools/

WORKDIR ${HOME}/tools/

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

## Install samtools
USER ubuntu
WORKDIR ${HOME}/tools/
RUN wget https://github.com/samtools/samtools/releases/download/1.2/samtools-1.2.tar.bz2 && \
    tar -xjf samtools-1.2.tar.bz2 && \
    rm -f samtools-1.2.tar.bz2 && \
    mv samtools-1.2 ${HOME}/tools/samtools

USER root
RUN cd ${HOME}/tools/samtools && \
    make -j && \
    make install

USER ubuntu

## Install fpfilter-tool
WORKDIR ${HOME}/tools/
RUN wget https://github.com/ucscCancer/fpfilter-tool/archive/master.zip && \
    unzip master.zip && \
    rm -f master.zip && \
    mv fpfilter-tool-master ${HOME}/tools/fpfilter-tool

WORKDIR ${HOME}
