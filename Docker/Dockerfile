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
RUN wget http://s3.amazonaws.com/plink2-assets/alpha2/plink2_linux_x86_64.zip && unzip plink2*.zip -d /usr/local/bin/ &&  rm plink2*.zip

ADD Scripts ../Scripts
ADD Data ../Data