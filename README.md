# dBg-min: linear-time minimization of the de Bruijn graph

Author: Nicola Prezza. 

### Description

Builds the de Bruijn graph from a fasta/fastq and computes in linear time the minimum equivalent Wheeler DFA (note: this may be larger than the minimum equivalent DFA). 

### Download

To clone the repository, run:

> git clone http://github.com/nicolaprezza/dBg-min

### Compile

The library has been tested under linux using g++ 10.3.0. 

We use cmake to generate the Makefile. Create a build folder in the main dBg-min folder:

> mkdir build

run cmake:

> cd build; cmake ..

and compile:

> make

### Run

After compiling, run 

>  dBg-min [-l nlines] [-a] input k

To build and minimize the de Bruijn graph of order k on the file "input" (a fastq file by default, or a fasta file if option -a is specified). if option -l nlines is specified, build the graph using only the first nlines sequences from the input file. 