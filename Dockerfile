FROM python:3.6-stretch

RUN apt-get update \
  && apt-get clean autoclean \
  && apt-get autoremove -y \
  && rm -rf /var/lib/apt/lists/*

# Copy source
COPY . variant-filtration-tool/
WORKDIR /variant-filtration-tool

# Install
RUN pip3 install .
