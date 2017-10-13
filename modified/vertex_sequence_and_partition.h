/*================================================================================
  
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

==================================================================================*/

#ifndef vertex_sequence_and_partition_h
#define vertex_sequence_and_partition_h

#include "weighted_graph.h"

typedef struct
{
    long long *sequence;
    long long number_of_subsets;
    long long *subset_size;
} sequence_and_partition;

sequence_and_partition * coloring_weighted(weighted_graph *graph,long long limit);
sequence_and_partition * coloring_unweighted(weighted_graph *graph,long long limit);

#endif
