FROM ubuntu:18.04
ENV TZ=Europe/Helsinki
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone
RUN apt-get update &&  apt-get install mercurial --yes && apt-get install zlib1g-dev libcurl4-openssl-dev libbz2-dev liblzma-dev --yes
RUN apt-get install wget make  gcc g++ datamash jq --yes unzip
RUN apt-get install python3 python3-numpy python3-pip  python3-pandas --yes
RUN pip3 install numpy --upgrade

# SAMTOOLS
RUN wget https://github.com/samtools/htslib/releases/download/1.10/htslib-1.10.tar.bz2 && \
tar -xvjf htslib-1.10.tar.bz2 && cd htslib-1.10 && ./configure && make && make install && cd .. && rm  htslib-1.*tar*

# BCFTOOLS
RUN wget https://github.com/samtools/bcftools/releases/download/1.10.2/bcftools-1.10.2.tar.bz2 && tar -xvjf bcftools-1.10.2.tar.bz2 && cd bcftools-1.10.2 && \
./configure && make && make install && cd ..  && rm  bcftools-1.10.2.tar.bz2 && export BCFTOOLS_PLUGINS=~/bcftools-1.10.2/plugins/

#PLINK
RUN wget https://s3.amazonaws.com/plink2-assets/plink2_linux_avx2_latest.zip && unzip plink2*.zip -d /usr/local/bin/ &&  rm plink2*.zip

RUN wget https://github.com/rgcgithub/regenie/releases/download/v2.0.2/regenie_v2.0.2.gz_x86_64_Linux_mkl.zip && unzip regen* && mv regenie_v2.0.2.gz_x86_64_Linux_mkl /usr/bin/regenie

ADD Scripts ../Scripts
ADD Data ../Data