ARG REGISTRY=docker.osdc.io/ncigdc
ARG BASE_CONTAINER_VERSION=latest

FROM ${REGISTRY}/python3.9-builder:${BASE_CONTAINER_VERSION} as builder

COPY ./ /gdc_filtration_tools

WORKDIR /gdc_filtration_tools

RUN pip install tox && tox -e build

RUN dnf install -y git

RUN git clone https://github.com/atks/vt.git <<EOF
cd vt
git submodule update --init --recursive
make
EOF

FROM ${REGISTRY}/python3.9-builder:${BASE_CONTAINER_VERSION}

LABEL org.opencontainers.image.title="gdc_filtration_tools" \
      org.opencontainers.image.description="This repository contains the source code used in the VCF variant filtration workflows within the GDC. A single CLI is generated with multiple subcommands." \
      org.opencontainers.image.source="https://github.com/NCI-GDC/variant-filtration-tool" \
      org.opencontainers.image.vendor="NCI GDC"

COPY --from=builder /gdc_filtration_tools/dist/*.whl /gdc_filtration_tools/
COPY requirements.txt /gdc_filtration_tools/

WORKDIR /gdc_filtration_tools

RUN pip install --no-deps --no-binary pysam -r requirements.txt \
	&& pip install --no-deps *.whl \
	&& rm -f *.whl requirements.txt

USER app

ENTRYPOINT ["gdc_filtration_tools"]

CMD ["--help"]
