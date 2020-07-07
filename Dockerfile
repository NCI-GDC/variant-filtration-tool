FROM python:3.6-stretch

RUN apt-get update \
  && apt-get clean autoclean \
  && apt-get autoremove -y \
  && rm -rf /var/lib/apt/lists/*

COPY ./dist /opt

RUN make -c /opt init-pip \
  && ln -s /opt/bin/gdc-filtration-tools /bin/gdc-filtration-tools

WORKDIR /opt

ENTRYPOINT ["/bin/gdc-filtration-tools"]

CMD ["--help"]
