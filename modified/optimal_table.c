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

static long long *mwc_solve_by_dp(weighted_graph *graph);

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
long long **create_optimal_table(long long *subset_size,long long number_of_subsets,weighted_graph *graph,long long limit)
{
    long long **optimal_table=(long long **)malloc(number_of_subsets * sizeof(long long *));

    for(long long i = 0; i < number_of_subsets; ++i) 
    {
        long long length=subset_size[i];
        long long *seq = (long long*)malloc(length*sizeof(long long));
        for(long long j = 0; j < length; ++j)
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
    graph: a vertex-weighted graph. the number of vertex must be less than sizeof(long long)
   <return>
    a part of the optimal table
 */
static long long *mwc_solve_by_dp(weighted_graph *graph)
{
    long long n = graph->n;
    long long *weight=graph->weight;
    long long **adjacency_matrix=get_bit_vector_adjacency_matrix(graph,graph->n);

    /*
       use [i][0] of adjacency_ matrix
       because graph->n is smaller than one word length
     */
    long long *adj0 = (long long*)malloc(n*sizeof(long long));
    adj0[0] = 0;
    for(long long j = 1; j < n; ++j)
    {
        adj0[j] = adjacency_matrix[j][0];
    }

    /* initialize table */
    long long *table = (long long*)malloc((1<<n)*sizeof(long long));
    table[0] = 0;

    /* dynamic programming */
    for(long long i = 0; i < n; ++i) 
    {
        long long start = 1 << i;
        long long end = 1 << (i+1);
        long long adji = adj0[i];
        long long weighti = weight[i];
        for(long long j = start; j < end; ++j) 
        {
            long long unused = table[j-start];
            long long used =  table[adji & j] + weighti;
            table[j] = unused < used ? used : unused;
        }
    }

    for(long long i = 0; i < n; ++i)
    {
        free(adjacency_matrix[i]);
    }
    free(adjacency_matrix);
    free(adj0);
    return table;
}
