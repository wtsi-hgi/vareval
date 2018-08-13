# build stage
FROM golang:alpine AS build-env
COPY . /src
RUN cd /src && go build -o mendErr

# final stage
#FROM vandhanak/bcftools:1.3.1 (works but old)
FROM ubuntu:xenial

# Needed by bcftools make install
RUN apt-get update &&  apt-get -y install build-essential && apt-get -y install zlib1g-dev && apt-get -y install  libbz2-dev \
 && apt-get -y install  liblzma-dev


COPY --from=build-env /src/bcftools-1.9 /tmp/bcftools-1.9
RUN cd /tmp/bcftools-1.9 \
  && ./configure --prefix=/usr/local \
  && make \
  && make install 

#Must be set to use plugins (that's changed location in 1.9)
ENV BCFTOOLS_PLUGINS=/usr/local/libexec/bcftools

COPY --from=build-env /src/mendErr /app/mendErr
WORKDIR /app

