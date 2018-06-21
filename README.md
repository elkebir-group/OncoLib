# OncoSim

OncoLib is a library for simulating the evolution and NGS sequencing of a metastatic tumor.

## Contents

  1. [Compilation instructions](#compilation)
     * [Dependencies](#dep)
     * [Compilation](#comp)
  2. [Usage instructions](#usage)
     * [I/O formats](#io)
       - [Clone tree](#clonetree)
       - [Vertex labeling](#vertexlabeling)
       - [Color map](#colormap)
       - [Frequency matrix](#frequencies)
       - [Read matrix](#reads)
     * [`simulate`](#simulate)
     * [Python module: `oncosim.so`](#oncosim.so)

<a name="compilation"></a>
## Compilation instructions

<a name="dep"></a>
### Dependencies

OncoSim is written in C++11 and thus requires a modern C++ compiler (GCC >= 4.8.1, or Clang). In addition, OncoSim has the following dependencies.

* [CMake](http://www.cmake.org/) (>= 3.0)
* [Boost](http://www.boost.org) (>= 1.55)
* [LEMON](http://lemon.cs.elte.hu/trac/lemon) graph library (>= 1.3)

[Graphviz](http://www.graphviz.org) is required to visualize the resulting DOT files, but is not required for compilation.

In case [doxygen](http://www.stack.nl/~dimitri/doxygen/) is available, extended source code documentation will be generated.

<a name="comp"></a>
### Compilation

To compile OncoSim, execute the following commands from the root of the repository:

    $ mkdir build
    $ cd build
    $ cmake ..
    $ make

In case CMake fails to detect LEMON, run the following command with adjusted paths:

    $ cmake -DLIBLEMON_ROOT=~/lemon ..

The compilation results in the following files in the `build` directory:

EXECUTABLE | DESCRIPTION
-----------|-------------
`simulate` | Simulate the migration and evolutionary history of a metastatic cancer.
`tree2dot` | Visualizes a given clone tree in DOT format.
`graph2dot` | Visualizes a given migration graph in DOT format.
`tree2freqs` | Extract mutation frequencies from a given clone tree.
`mix` | Generate mixtures of the leaves of a given clone tree.
`sequence`  | Generate a read matrix (containing NGS read counts) from a given clone tree.
`reads2freqs` | Converts reads into mutation frequency confidence intervals.
`downsample` | Down sample a given read matrix.
`oncosim.so` | Python library.

<a name="usage"></a>
## Usage instructions

<a name="io"></a>
### I/O formats

Below we describe the various formats used by OncoSim.

<a name="clonetree"></a>
#### Clone tree

The first line of a clone tree file indicates the number `t` of edges. The subsequent `t` lines each encode an edge, by listing the labels of the incident vertices separated by a space or tab character. Next, the number `s` of leaves is specified. The subsequent `s` lines each encode a leaf. Each line is a string whose values are separated by spaces or tabs. The first value is the leaf label, followed by the location. Next, the number of samples is specified. For each sample the mixture proportion is specified. The remaining values correspond to mutations.

    8 #edges
    26;28 26;28_M1
    25 26;28
    25 25_M1
    24 24_P
    0;2;3;4;5;6;7;8;9;10;11;12;13;14;15;16;17;18;19;20;21;22;23 25
    0;2;3;4;5;6;7;8;9;10;11;12;13;14;15;16;17;18;19;20;21;22;23 24
    0;2;3;4;5;6;7;8;9;10;11;12;13;14;15;16;17;18;19;20;21;22;23 22;23_P
    GL 0;2;3;4;5;6;7;8;9;10;11;12;13;14;15;16;17;18;19;20;21;22;23
    4 #leaves
    22;23_P P 2 1 0.850564 0 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23
    24_P P 2 0 0.149436 0 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24
    26;28_M1 M1 2 0 0.46333 0 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 25 26 28
    25_M1 M1 2 1 0.53667 0 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 25

<a name="vertexlabeling"></a>
#### Vertex labeling

A vertex labeling assigns a location to each vertex of a clone tree. Each line contains two values, the vertex label and the location label separated by a space or tab character. For example:

    GL P
    0 P
    1 P
    1_P P
    2;3;4;5;6;7;8;9;10;11;12;13;14;15;16;17;18;19;20;21;22;23 P
    22;23_P P
    24 P
    24_P P
    25_M1 M1
    26;28 M1
    26;28_M1 M1
    25 M1

<a name="colormap"></a>
#### Color map

A color map file assigns a color to each location. Each line contains two values, the location label and an integer indicating the color.

    P 1
    M_1 2

<a name="frequencies"></a>
#### Frequency matrix

A frequency file encodes the frequency of every mutation in a sample. It is a tab separated file. The first line lists the number of anatomical sites followed by the number of samples and then the number of mutations, each on separate lines. The fourth line is ignored but describes the format of the rest of the file. Each subsequent line encodes the cell frequency of a mutation in a sample: first the sample 0-based index is given, followed by the label of the sample, the 0-based index of the anatomical site, the anatomical site label, the 0-based index of the mutation, the label of the mutation, the frequency lower bound and upper bound.

    2 #anatomical sites
    4 #samples
    27 #characters
    #sample_index	sample_label	anatomical_site_index	anatomical_site_label	character_index	character_label	f-	f+
    0	M1_0	0	M1	0	0	0.962918	1
    0	M1_0	0	M1	1	2	0.892128	1
    0	M1_0	0	M1	2	3	0.870492	1
    0	M1_0	0	M1	3	4	0.895341	1
    0	M1_0	0	M1	4	5	0.775452	1
    ...

<a name="reads"></a>
#### Read matrix
A reads file encodes the number of variant and reference reads of every mutation in each sample. It is a tab separated file. The first line lists the number of anatomical sites followed by the number of samples and then the number of mutations, each on separate lines. The fourth line is ignored but describes the format of the rest of the file. Each subsequent line encodes the variant and reference read counts of a mutation in a sample: first the sample 0-based index is given, followed by the label of the sample, the 0-based index of the anatomical site, the anatomical site label, the 0-based index of the mutation, the label of the mutation, the number of reference reands and the number of variant reads.

    2 #anatomical sites
    4 #samples
    27 #characters
    #sample_index	sample_label	anatomical_site_index	anatomical_site_label	character_index	character_label	ref	var
    0	M1_0	0	M1	0	0	72	91
    0	M1_0	0	M1	1	2	83	90
    0	M1_0	0	M1	2	3	77	81
    0	M1_0	0	M1	3	4	75	83
    0	M1_0	0	M1	4	5	83	72
    ...

<a name="simulate"></a>
### `simulate`

    $ ../build/simulate -h
    Usage:
      ../build/simulate [--help|-h|-help] [-C int] [-D num] [-E num] [-K num]
         [-N int] [-P num] [-c str] [-f num] [-k int] [-kP int] [-m int]
         [-mig num] [-mut num] [-o str] [-p int] [-s int] [-v]
    Where:
      --help|-h|-help
         Print a short help message
      -C int
         Target coverage (default: 200)
      -D num
         Driver probability (default: 1e-7)
      -E num
         Per base sequencing error rate (default: 0)
      -K num
         Carrying capacity (default: 5e4)
      -N int
         Number of successful simulations (default: -1)
      -P num
         Purity (default: 1)
      -c str
         Color map file
      -f num
         Mutation frequency threshold (default: 0.05)
      -k int
         Number of samples per anatomical site (default: 2)
      -kP int
         Number of samples for the primary tumor (default: 2)
      -m int
         Maximum number of detectable anatomical sites (default: 8)
      -mig num
         Migration rate (default: 1e-6)
      -mut num
         Mutation rate (default: 0.1)
      -o str
         Output directory
      -p int
         Allowed migration patterns:
           0 : mS (default)
           1 : mS, S
           2 : mS, S, M
           3 : mS, S, M, R
      -s int
         Random number generator seed (default: 0)
      -v
         Verbose


<a name="oncosim.so"></a>
### Python module: `oncosim.so`

The Python module `oncosim.so` has the following eight functions:

1. simulate
2. tree2dot
3. graph2dot
4. tree2freqs
5. mix
6. sequence
7. reads2freqs
8. downsample

These functions correspond to the standalone executables described previously. Please refer to [example/example.ipynb](example/example.ipynb) for example usage.
