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
    sqlite3 \
    git-core

RUN apt-get clean \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

ENV variant-filtration-tool 2.6

WORKDIR /opt

## Install vt
## Only needed for PINDEL
#RUN git clone https://github.com/atks/vt.git && \
#    cd vt && \
#    git checkout -b 40d79cfb6e0ccbfca93436afb81b2abbc7aa1c26 && \
#    make

## Install variant-filtration-tool
RUN mkdir /opt/variant-filtration-tool
WORKDIR /opt/variant-filtration-tool
ADD utils /opt/variant-filtration-tool 
ADD requirements.txt /opt/variant-filtration-tool 

RUN cd /opt/variant-filtration-tool && \
    pip3 install -r ./requirements.txt

WORKDIR /opt 
