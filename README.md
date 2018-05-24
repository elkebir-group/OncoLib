# OncoSim

OncoSim is a library for simulating the evolution and NGS sequencing of a metastatic tumor.

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

<a name="frequencies"></a>
#### Frequency matrix

<a name="reads"></a>
#### Read matrix
A frequency file encodes the frequency of every mutation (cluster) in an anatomical site (sample). It is a tab separated file. The first line lists the number of anatomical sites followed by the number of samples and then the number of mutations, each on separate lines. The fourth line is ignored but describes the format of the rest of the file. Each subsequent line encodes the cell frequency of a mutation in a sample: first the sample 0-based index is given, followed by the label of the sample, the 0-based index of the anatomical site, the anatomical site label, the 0-based index of the mutation, the label of the mutation, the frequency lower bound and upper bound.

    6 #anatomical sites							
    6 #samples							
    10 #mutation clusters							
    #sample_index	sample_label	anatomical_site_index	anatomical_site_label	character_index	character_label	f_lb	f_ub
    0	breast	0	breast	0	1	0.503628522	0.545237495
    0	breast	0	breast	1	2	0	0.01213794

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

Please refer to [example/example.ipynb](example/example.ipynb) for example usage.
