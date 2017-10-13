The original code by Shimizu and Yamaguchi is in the `original-code` directory.
This was downloaded from http://www2.kobe-u.ac.jp/%7Eky/otclique/otclique.html
on 13 October 2017.

This modified version uses `long long` rather than `int`, to enable the use
of larger weights.  On my machine, this makes the code something like 30% slower
than the original program.  Note that the change of data types was made using
a simple text substitution; see the `scripts` directory.

The modified program also allows a time limit in seconds to be set using a
third command-line argument.  To use the default value of the l parameter along
with a timeout, use -1 for the second command-line argument.

This code is similar, but not identical, to the code used for the CP 2017
paper *On Maximum Weight Clique Algorithms, and How They Are Evaluated*.  The
used for that paper was based on an earlier version of OTClique provided
by the program's authors.
