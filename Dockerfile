FROM python:3.6-stretch

RUN apt-get update \
  && apt-get clean autoclean \
  && apt-get autoremove -y \
  && rm -rf /var/lib/apt/lists/*

COPY ./dist /opt

RUN make -C /opt init-pip \
  && ln -s /opt/bin/gdc_filtration_tools /bin/gdc_filtration_tools

WORKDIR /opt

ENTRYPOINT ["/bin/gdc-filtration-tools"]

CMD ["--help"]
