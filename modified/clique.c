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

#include "clique.h"
#include "weighted_graph.h"
#include <assert.h>

long long is_clique(clique *clq,weighted_graph *graph)
{
    long long **adj=graph->adjacency_matrix;
    long long size=clq->size;
    long long *set=clq->set;
    for(long long i=0;i<size-1;i++)
    {
        for(long long j=i+1;j<size;j++)
        {
            if(!adj[set[i]][set[j]])
            {
                return 0;
            }
        }
    }
    return 1;
}
