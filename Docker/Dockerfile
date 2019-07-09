FROM ubuntu:18.04
ENV TZ=Europe/Helsinki
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone
RUN apt-get update
RUN apt-get install mercurial --yes && apt-get install gcc zlib1g-dev wget make libcurl4-openssl-dev libbz2-dev liblzma-dev --yes
RUN apt-get install bzip2 curl unzip gawk python3 make gcc g++ --yes
RUN apt-get install libz-dev libbz2-dev liblzma-dev littler libssl-dev python3-numpy python3-pandas python3-pip python3-scipy --yes
RUN pip3 install scipy && pip3 install numpy --upgrade

# SAMTOOLS
RUN wget https://github.com/samtools/htslib/releases/download/1.9/htslib-1.9.tar.bz2 && \
tar -xvjf htslib-1.9.tar.bz2 && cd htslib-1.9 && ./configure && make && make install && cd .. && rm  htslib-1.9.tar.bz2

RUN wget https://github.com/samtools/bcftools/releases/download/1.9/bcftools-1.9.tar.bz2 && tar -xvjf bcftools-1.9.tar.bz2 && cd bcftools-1.9 && \
./configure && make && make install && cd ..  && rm  bcftools-1.9.tar.bz2 && export BCFTOOLS_PLUGINS=~/bcftools-1.9/plugins/


ADD Scripts ../Scripts