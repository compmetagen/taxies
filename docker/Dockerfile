FROM debian:9
MAINTAINER Davide Albanese <davide.albanese@gmail.com>

RUN apt-get update && apt-get install -y \
    build-essential \
    pkg-config \
    automake \
    wget \
    git \
    zip \
    tar \
    gfortran \
    libblas-dev \
    liblapack-dev \
    libfreetype6 \
    libfreetype6-dev \
    libpng-dev \
    liblzma-dev \
    libbz2-dev \
    zlib1g \
    zlib1g-dev \
    python3 \
    python3-dev \
    python3-pip \
    python3-numpy \
    python3-scipy \
    python3-matplotlib \
  && rm -rf /var/lib/apt/lists/*

RUN pip3 install 'Click==6.7' 'pandas==0.19.2' 'dendropy==4.2.0' \
    'biopython==1.68'

# Prodigal v2.6.3
RUN wget https://github.com/hyattpd/Prodigal/archive/v2.6.3.tar.gz \
    -q -O /tmp/v2.6.3.tar.gz \
  && tar -zxf /tmp/v2.6.3.tar.gz -C /tmp
WORKDIR /tmp/Prodigal-2.6.3
RUN make install
WORKDIR /
RUN rm -rf /tmp/v2.6.3.tar.gz /tmp/Prodigal-2.6.3

# HMMER v3.2.1
RUN wget http://eddylab.org/software/hmmer/hmmer-3.2.1.tar.gz \
    -q -O /tmp/hmmer-3.2.1.tar.gz \
  && tar -zxf /tmp/hmmer-3.2.1.tar.gz -C /tmp
WORKDIR /tmp/hmmer-3.2.1
RUN ./configure && make && make install
WORKDIR /tmp/hmmer-3.2.1/easel
RUN make install
WORKDIR /
RUN rm -rf /tmp/hmmer-3.2.1.tar.gz /tmp/hmmer-3.2.1

# VSEARCH v2.9.1
RUN wget https://github.com/torognes/vsearch/archive/v2.9.1.tar.gz \
    -q -O /tmp/v2.9.1.tar.gz \
  && tar -zxf /tmp/v2.9.1.tar.gz -C /tmp
WORKDIR /tmp/vsearch-2.9.1
RUN ./autogen.sh && ./configure && make && make install
WORKDIR /
RUN rm -rf /tmp/v2.9.1.tar.gz /tmp/vsearch-2.9.1

# Taxies
RUN git clone https://github.com/compmetagen/taxies.git /tmp/taxies/
WORKDIR /tmp/taxies/
RUN python3 setup.py install
WORKDIR /
RUN rm -rf /tmp/taxies/
