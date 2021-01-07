# DNA_stability

[![Build Status](https://travis-ci.com/jamesmtuck/DNA_stability.svg?token=rCvdBqMzwWyNvxxUUbSh&branch=main)](https://travis-ci.com/jamesmtuck/DNA_stability)

- [Overview](#overview)
- [Documentation](#documentation)
- [System Requirements](#system-requirements)
- [Installation Guide](#installation-guide)
- [License](#license)
- [Issues](https://github.com/jamesmtuck/DNA_stability/issues)

# Overview

Analysis of DNA stability on design of error correction codes for DNA-based data storage. 

# Documentation

As documentation for the softwarwe becomes available, it will be placed under the docs folder.

# System Requirements

## Hardware Requirements
This code requires only a standard computer with enough RAM and compute power to support the needed operations.

## Software Requirements
### OS Requirements
This package is supported for macOS and Linux. The package has been tested on the following systems:

+ macOS: Catalina 10.15.3
+ Linux: Ubuntu 18.04.3

Note that most OSes will support our software by using Docker.

### Other Dependences

The python packages used for high precision floating point operations have dependences with libraries written in other languages. They can be installed using common package managers:
```
sudo apt-get install -qq libgmp-dev libmpc-dev libmpfr-dev
```
### Python Dependences

Our code has been tested on python versions 3.5 to 3.8. It has the following dependences:

```
nose
sphinx
biopython
editdistance
statistics
matplotlib
numpy
scipy
gmpy2
```

# Installation Guide

## Use your local python environment
If you already have python 3 installed on your system, the simplest thing to do is download or checkout the code from GitHub.  Then in the DNA_stability directory, run the following commands:

    git clone https://github.com/jamesmtuck/DNA_stability
    cd DNA_stability
    pip install -r requirements.txt

## Use Docker

If you do not already have Docker, you will need to install Docker on your system. It is available for free for most versions of Windows, Linux, and MacOS. You may need to be the owner or administrator of the system to install Docker.

Instructions for setting up Docker.  From a command prompt, run these commands:

    git clone https://github.com/jamesmtuck/DNA_stability
    cd DNA_stability
    docker build -t stability:1.0 .
    docker run -it -v `pwd`:/stability stability:1.0 /bin/bash

This will bring up a command prompt in a Linux container where commands can be executed. 

## Run Our Analyses

Run the main script this way:

    python3 stability.py
    
This command should run to completion within a few minutes.

# License

This software is released under the MIT License.


