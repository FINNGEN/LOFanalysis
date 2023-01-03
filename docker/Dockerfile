FROM eu.gcr.io/finngen-refinery-dev/bioinformatics:0.8

RUN apt-get update && apt-get upgrade -y

ADD ./docker/requirements.txt ../requirements.txt
RUN pip3 install -r requirements.txt

ADD scripts ../scripts
ADD data ../data