FROM ubuntu:18.04
ENV TZ=Europe/Helsinki
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone
RUN apt-get update &&  apt-get install mercurial --yes && apt-get install zlib1g-dev libcurl4-openssl-dev libbz2-dev liblzma-dev --yes
RUN apt-get install wget make  gcc g++ datamash --yes
RUN apt-get install python3 python3-numpy python3-pip  python3-pandas --yes
RUN pip3 install numpy --upgrade


# SAMTOOLS
RUN wget https://github.com/samtools/htslib/releases/download/1.9/htslib-1.9.tar.bz2 && \
tar -xvjf htslib-1.9.tar.bz2 && cd htslib-1.9 && ./configure && make && make install && cd .. && rm  htslib-1.9.tar.bz2

RUN wget https://github.com/samtools/bcftools/releases/download/1.9/bcftools-1.9.tar.bz2 && tar -xvjf bcftools-1.9.tar.bz2 && cd bcftools-1.9 && \
./configure && make && make install && cd ..  && rm  bcftools-1.9.tar.bz2 && export BCFTOOLS_PLUGINS=~/bcftools-1.9/plugins/


ADD Scripts ../Scripts