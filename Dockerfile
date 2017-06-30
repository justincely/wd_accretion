# Base image to use
FROM ubuntu:xenial

MAINTAINER Justin Ely <justincely@gmail.com>

RUN apt-get update
RUN apt-get upgrade -y
RUN apt-get install -y git wget bzip2 gfortran

# Get the mercury code and simulation support code
COPY . /accretion_sim

# get conda
RUN wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
RUN bash miniconda.sh -b -p /root/miniconda
ENV PATH="/root/miniconda/bin:${PATH}"

RUN conda install numpy matplotlib astropy

# 
WORKDIR accretion_sim/mercury_code
RUN gfortran -o mercury6 mercury6_3.for
RUN gfortran -o element6 element6.for
RUN gfortran -o close6 close6.for

ENTRYPOINT ["./generate_system"]