FROM python:3.6

# Copy source
COPY ./utils variant-filtration-tool/
ADD requirements.txt variant-filtration-tool/
WORKDIR /variant-filtration-tool

# Install
RUN pip3 install -r requirements.txt
