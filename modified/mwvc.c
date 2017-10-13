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
    long long limit;
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

    weighted_graph *complement_graph = get_complement_graph(input_graph);
    clique *maximum_weight_clique=otclique(complement_graph,limit);

    long long all_weight_sum = 0;
    {
        long long *weight=input_graph->weight;
        for(long long i = 0; i<input_graph->n; i++)
        {
            all_weight_sum += weight[i];
        }
    }

    long long *mwc = (long long *)calloc(input_graph->n,sizeof(long long));

    for(long long i = 0; i < maximum_weight_clique->size; ++i)
    {
        mwc[ maximum_weight_clique->set[i] ] = 1;
    }

    printf("Minimum weight = %lld\n", all_weight_sum-maximum_weight_clique->weight);
    printf("The minimum weight vertex cover has %lld vertices,\n [",input_graph->n - maximum_weight_clique->size);
    for(long long i = 0; i < input_graph->n; ++i)
    {
        if(mwc[i] == 0)
        {
            printf(" %lld", i+1);
        }
    }
    printf(" ]\n");

    assert(is_clique(maximum_weight_clique,complement_graph));

    free(input_graph->adjacency_matrix);
    free(input_graph->weight);
    free(input_graph);
    free(complement_graph->adjacency_matrix);
    free(complement_graph->weight);
    free(complement_graph);
    free(maximum_weight_clique->set);
    free(maximum_weight_clique);
    free(mwc);
    return 0;
}
