FROM python:3.6

# Copy source
COPY . variant-filtration-tool/
WORKDIR /variant-filtration-tool

# Install
RUN pip3 install .
