
This program is an implementation of the algorithm described in the paper "Directed, weighted and overlapping benchmark graphs for community detection algorithms", written by Andrea Lancichinetti and Santo Fortunato. In particular, this program is to produce binary networks with overlapping nodes. 
Each feedback is very welcome. If you have found a bug or have problems, or want to give advises, please contact us:

andrea.lancichinetti@isi.it
fortunato@isi.it


Turin, 29 October 2009
---------------------------------------------------------------





-------------------- How to compile -----------------------------------
In order to compile, type:

make

-------------------- How to run the program ---------------------------


To run the program, type:

./benchmark [FLAG] [P]


[FLAG]		[P]

-N		number of nodes
-k		average degree
-maxk		maximum degree
-mu		mixing parameter
-t1		minus exponent for the degree sequence
-t2		minus exponent for the community size distribution
-minc		minimum for the community sizes
-maxc		maximum for the community sizes
-on		number of overlapping nodes
-om		number of memberships of the overlapping nodes
-C              [average clustering coefficient]
----------------------

In this program you can assign the number of overlapping nodes (option -on) and assign the number of memberships for them (option -om). The other nodes will have only one membership. For instance, typing 

./benchmark [flags...] -N 1000 -on 20 -om 2

will produce a network with 1000 nodes, 980 with only one membership and 20 nodes with two memberships.

It is also possible to set the parameters writing flags and relative numbers in a file (look at flags.dat to see an example). To specify the file, use the option:

-f	[filename]

You can set the parameters both writing some of them in the file, and using flags from the command line for others (if you set a parameter twice, the latter one will be taken).

-N, -k, -maxk, -mu have to be specified. For the others, the program can use default values:

t1=2, 	t2=1, 	on=0,	om=0,	minc and maxc will be chosen close to the degree sequence extremes.

----------------------

The flag -C is not mandatory. If you use it, the program will perform a number of rewiring to increase the average cluster coefficient up to the wished value.
Since other constraint must be fulfilled, if the wished value will not be reached after a certain time, the program will stop (displaying a warning).




-------------------- Other options ---------------------------

To have a random network use:
-rand

Using this option will set mu=0, and minc=maxc=N, i.e. there will be one only community.

Use option:
-sup (-inf) 

if you want to produce a benchmark whose distribution of the ratio of external degree/total degree is superiorly (inferiorly) bounded by the mixing parameter. In other words, if you use one of these options, the mixing parameter is not the average ratio of external degree/total degree (as it used to be) but the maximum (or the minimum) of that distribution. When using one of these options, what the program essentially does is to approximate the external degree always by excess (or by defect) and if necessary to modify the degree distribution. Nevertheless, this last possibility occurs for a few nodes and numerical simulations show that it does not affect the degree distribution appreciably.


-------------------- Examples ---------------------------

Example1:
./benchmark -N 1000 -k 15 -maxk 50 -mu 0.1 -minc 20 -maxc 50
Example2:
./benchmark -f flags.dat -t1 3

If you want to produce a kind of Girvan-Newman benchmark, you can type:
./benchmark -N 128 -k 16 -maxk 16 -mu 0.1 -minc 32 -maxc 32


-------------------- Output ---------------------------

Please note that the community size distribution can be modified by the program to satisfy several constraints (a warning will be displayed).

The program will produce three files:

1) network.dat contains the list of edges (nodes are labelled from 1 to the number of nodes; the edges are ordered and repeated twice, i.e. source-target and target-source).

2) community.dat contains a list of the nodes and their membership (memberships are labelled by integer numbers >=1).

3) statistics.dat contains the degree distribution (in logarithmic bins), the community size distribution, and the distribution of the mixing parameter.



-------------------- Seed for the random number generator ---------------------------
 

-In the file time_seed.dat you can edit the seed which generates the random numbers. After reading, the program will increase this number by 1 (this is done to generate different networks running the program again and again). If the file is erased, it will be produced by the program again.


