Copyright (c) 2017 Satoshi SHIMIZU and Kazuaki YAMAGUCHI. All rights reserved.

< Academic use >
    1. Any modification is allowed for academic purposes.
       If you make a modified version, please redistribute its source code.

    2. Please cite our paper when using this code.

       @article{otclique,
                title={Fast maximum weight clique extraction algorithm: Optimal tables for branch-and-bound},
                author={Shimizu, Satoshi and Yamaguchi, Kazuaki and Saitoh, Toshiki and Masuda, Sumio},
                journal={Discrete Applied Mathematics},
                volume={223},
                pages={120--134},
                year={2017},
                publisher={Elsevier}
       }

< Commercial use >
    Basically, this code can be used only for academic purposes.
    For commercial use, please make a contact with us.

< Usage >
    Run command as following.

    otclique input_file [subset_size_limit]

    Input file must be an ascii file of the DIMACS format.
    See [ ftp://dimacs.rutgers.edu/pub/challenge/graph/doc/ccformat.dvi ].
    In addition, lines of "n V WEIGHT" denotes the vertex V has weight W.
    V and W must be positive integer values.
    The default value of the weight is 1.

    The optional argument subset_size_limit is a parameter used in the algorithm.
    It must be less than 32.  Set it to appropriate value for the amount of memory.
    If number of vertices is less than 1500, default value is 20, otherwise it is set to 25 by default.
