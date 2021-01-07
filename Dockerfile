# Get repo with latest python repos
FROM python:latest

LABEL maintainer="jtuck@ncsu.edu"

WORKDIR /DNA_stability

COPY . .

# Install libs
RUN apt-get update       \
  &&  apt-get install -y \
     libgmp-dev   \
     libmpc-dev          \
     libmpfr-dev         \
  && apt-get clean
    

# Install the needed tools
RUN pip3 --no-cache-dir install -r requirements.txt

# Set python path to include the code in the doris module
ENV "PYTHONPATH" "/DNA_stability/"
