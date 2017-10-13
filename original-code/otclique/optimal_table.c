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

#include <stdlib.h>
#include "weighted_graph.h"

static int *mwc_solve_by_dp(weighted_graph *graph);

/*
   create optimal table
   <args>
    subset_size: the size of each subset.
    number_of_subsets: the number of subsets.
    graph: a graph. vertices must be already reordered and re-indexed.
    limit: the limit size of each subsets
   <return>
    the optimal table
 */
int **create_optimal_table(int *subset_size,int number_of_subsets,weighted_graph *graph,int limit)
{
    int **optimal_table=(int **)malloc(number_of_subsets * sizeof(int *));

    for(int i = 0; i < number_of_subsets; ++i) 
    {
        int length=subset_size[i];
        int *seq = (int*)malloc(length*sizeof(int));
        for(int j = 0; j < length; ++j)
        {
            seq[j] = limit * i + j;
        }
        weighted_graph *graph2 = create_vertex_induced_subgraph(seq,length,graph);
        optimal_table[i]=mwc_solve_by_dp(graph2);

        free(graph2->adjacency_matrix);
        free(graph2->weight);
        free(graph2);
        free(seq);
    }
    return optimal_table;
}

/*
   calculate all exact solutions of all subgraphs of givin graph
   <args>
    graph: a vertex-weighted graph. the number of vertex must be less than sizeof(int)
   <return>
    a part of the optimal table
 */
static int *mwc_solve_by_dp(weighted_graph *graph)
{
    int n = graph->n;
    int *weight=graph->weight;
    int **adjacency_matrix=get_bit_vector_adjacency_matrix(graph,graph->n);

    /*
       use [i][0] of adjacency_ matrix
       because graph->n is smaller than one word length
     */
    int *adj0 = (int*)malloc(n*sizeof(int));
    adj0[0] = 0;
    for(int j = 1; j < n; ++j)
    {
        adj0[j] = adjacency_matrix[j][0];
    }

    /* initialize table */
    int *table = (int*)malloc((1<<n)*sizeof(int));
    table[0] = 0;

    /* dynamic programming */
    for(int i = 0; i < n; ++i) 
    {
        int start = 1 << i;
        int end = 1 << (i+1);
        int adji = adj0[i];
        int weighti = weight[i];
        for(int j = start; j < end; ++j) 
        {
            int unused = table[j-start];
            int used =  table[adji & j] + weighti;
            table[j] = unused < used ? used : unused;
        }
    }

    for(int i = 0; i < n; ++i)
    {
        free(adjacency_matrix[i]);
    }
    free(adjacency_matrix);
    free(adj0);
    return table;
}
