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

#include "weighted_graph.h"
#include "clique.h"
#include "otclique.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

int main(int argc, char *argv[])
{
    int limit;
    weighted_graph *input_graph;
    switch (argc)
    {
        case 3 :
            input_graph=read_graph(argv[1]);
            limit=atoi(argv[2]);
            break;
        case 2 :
            input_graph=read_graph(argv[1]);
            if( input_graph->n <=1500 )
            {
                limit=25;
            }
            else
            {
                limit=20;
            }
            break;
        default:
            fprintf(stderr,"Usage: %s file [subset_size_limit]\n",argv[0]);
            return 1;
    }

    clique *maximum_weight_clique=otclique(input_graph,limit);

    printf("Maximum weight = %d\n", maximum_weight_clique->weight);
    printf("The maximum weight clique has %d vertices,\n [",maximum_weight_clique->size);
    for(int i = 0; i < maximum_weight_clique->size; ++i)
    {
        printf(" %d", maximum_weight_clique->set[i]+1);
    }
    printf(" ]\n");

    assert(is_clique(maximum_weight_clique,input_graph));

    free(input_graph->adjacency_matrix);
    free(input_graph->weight);
    free(input_graph);
    free(maximum_weight_clique->set);
    free(maximum_weight_clique);
    return 0;
}
