FROM ubuntu:16.04

RUN apt-get update && apt-get install mercurial --yes && \
apt-get install gcc zlib1g-dev wget make libcurl4-openssl-dev libbz2-dev liblzma-dev \
bzip2 curl unzip gawk python python3 make gcc g++ libz-dev libbz2-dev liblzma-dev libreadline6 libreadline6-dev littler libssl-dev python3-numpy python3-pandas python3-pip --yes &&\
pip3 install networkx && pip3 install numpy --upgrade


# PLINK 1.9
RUN wget http://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20181202.zip && unzip plink_linux_x86_64_20181202.zip -d /usr/local/bin/ &&  rm plink*.zip

#PLINK2
RUN wget http://s3.amazonaws.com/plink2-assets/alpha1/plink2_linux_x86_64.zip  && unzip plink2*.zip -d /usr/local/bin/ &&  rm plink2*.zip

ADD Scripts ../Scripts
