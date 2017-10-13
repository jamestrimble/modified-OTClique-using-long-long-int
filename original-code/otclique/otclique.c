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
#include "otclique.h"
#include "weighted_graph.h"
#include "optimal_table.h"
#include "vertex_sequence_and_partition.h"
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <memory.h>
#include <limits.h>

static clock_t start_precomputation;
static clock_t end_precomuputation;
static clock_t start_branch_and_bound;
static clock_t end_branch_and_bound;
static unsigned long branch_count;
static int **adjacency_matrix;
static int *weight;
static int *c;
static int *record;
static int record_weight;
static int record_size;
static int *current;
static int current_size;
static int current_weight;
static int **optimal_table;
static int limit;
static int *msb_table;
static int weighted;
static int number_of_subsets,*subset_size;
static weighted_graph *input_graph;
static weighted_graph *reconstructed_graph;
static int *seq=NULL;

static void precomputation();
static void branch_and_bound();
static void expand(int *set,int set_size,int upper);

clique * otclique(weighted_graph *graph,int subset_size_limit)
{
    double sec_precomputation;
    double sec_branch_and_bound;
    double sec_total;
    clique *maximum_weight_clique;
    int n=graph->n;
    branch_count=0;
    input_graph=graph;
    limit=subset_size_limit;

    printf("Subset size limit = %d\n", limit);
    start_precomputation=clock();

    /* check if the graph is weighted or unweighted */
    {
        weighted=0;
        int weight0=graph->weight[0];
        int *weight=graph->weight;
        for(int i=0;i<n;i++)
        {
            if(weight0 != weight[i])
            {
                weighted=1;
                break;
            }
        }
    }

    /* create msb_table to get msb of any bit vector */
    msb_table = (int *)malloc((1 << limit)*sizeof(int));
    msb_table[0] = -1;
    for(int i = 0; i < limit; ++i) 
    {
        int from = 1 << i;
        int to = 1 << (i + 1);
        for(int j = from; j < to; ++j) 
        {
            msb_table[j] = i;
        }
    }

    /* precomputation phase */
    precomputation();

    end_precomuputation=clock();
    sec_precomputation=(double)(end_precomuputation-start_precomputation)/CLOCKS_PER_SEC;
    printf("%d subsets created from %d vertices \n", number_of_subsets, n);
    /* print record */
    printf("Precomputation phase = %.2f sec.\n",
            sec_precomputation);

    start_branch_and_bound=clock();
    /* branch-and-bound phase */
    branch_and_bound();

    end_branch_and_bound=clock();
    sec_branch_and_bound=(double)(end_branch_and_bound-start_branch_and_bound)/CLOCKS_PER_SEC;
    sec_total=(double)(end_branch_and_bound-start_precomputation)/CLOCKS_PER_SEC;

    /* print record */
    printf("Branch-and-bound phase = %.2f sec.\n",
            sec_branch_and_bound);
    printf("Branch-and-bound iterations = %ld (recursive calls)\n",
            branch_count);
    printf("Total time = %.2f sec.\n", sec_total);

    maximum_weight_clique=(clique *)malloc(sizeof(clique));
    maximum_weight_clique->size=record_size;
    maximum_weight_clique->weight=record_weight;
    maximum_weight_clique->set=(int *)malloc(sizeof(int) * record_size);
    for(int i = 0; i < record_size; ++i)
    {
        maximum_weight_clique->set[i]=seq[record[i]];
    }

    for(int i = 0;i<number_of_subsets;i++)
    {
        free(optimal_table[i]);
    }
    free(optimal_table);
    free(msb_table);
    {
        int n_r=reconstructed_graph->n;
        for(int i=0; i<n_r; i++)
        {
            free(adjacency_matrix[i]);
        }
    }
    free(adjacency_matrix);
    free(reconstructed_graph->adjacency_matrix);
    free(reconstructed_graph);
    free(seq);
    free(subset_size);
    free(record);
    free(weight);
    return maximum_weight_clique;
}

/*
   Precomputation phase.
   1. Make a vertex sequence and partition.
   2. Make the optimal tables.
 */
static void precomputation()
{
    /* make a vertex sequence and partition */
    sequence_and_partition *seq_and_partition;
    if(weighted)
    {
        seq_and_partition=coloring_weighted(input_graph,limit);
    }
    else
    {
        seq_and_partition=coloring_unweighted(input_graph,limit);
    }

    int *seq0=seq_and_partition->sequence;
    number_of_subsets=seq_and_partition->number_of_subsets;
    subset_size = seq_and_partition->subset_size;
    free(seq_and_partition);

    /* encode to bit vector*/
    seq=(int *)calloc(number_of_subsets*limit,sizeof(int));
    {
        int k=0;
        for(int i=0;i<number_of_subsets;i++)
        {
            int h=i*limit;
            for(int l=0;l<subset_size[i];l++)
            {
                seq[h++]=seq0[k++];
            }
        }
    }
    free(seq0);
    reconstructed_graph=create_vertex_induced_subgraph(seq,number_of_subsets*limit,input_graph);
    adjacency_matrix=get_bit_vector_adjacency_matrix(reconstructed_graph,limit);
    weight=reconstructed_graph->weight;

    /* create optimal tables */
    optimal_table=create_optimal_table(subset_size,number_of_subsets,reconstructed_graph,limit);
}

/*
   Branch-and-bound phase.
 */
static void branch_and_bound()
{
    int n=input_graph->n;

    /* initialize variables used in branch-and-bound phase */
    record=(int *)malloc(sizeof(int)*n);
    record_weight=0;
    record_size=0;
    current = (int *)malloc(sizeof(int)*n);
    c = (int *)malloc(sizeof(int) * (number_of_subsets*limit));
    for(int i=0;i<number_of_subsets*limit;i++)
    {
        c[i]=INT_MAX/2;
    }

    int stop=n;
    if(weighted)
    {
        stop=n*0.8;
    }

    /* main loop */
    int *set = (int *)calloc(number_of_subsets,sizeof(int));
    {
        int i=0;int j=0;int l=0;
        for(i=0; i<number_of_subsets; i++)
        {
            for(;j<subset_size[i];j++)
            {
                if(l++==stop)
                {
                    goto nobs; /* stop calculation of c[] */
                }
                set[i] += (1<<j);
                int v = (i*limit) + j;
                int* adjv=adjacency_matrix[v];
                int* set2=(int *)calloc(i+1,sizeof(int));
                int upper=0;
                int k=(v-1)/limit+1;
                while(k--)
                {
                    set2[k] = set[k] & adjv[k];
                    upper+=optimal_table[k][set2[k]];
                }
                current_size=1;
                current_weight=weight[v];
                current[0]=v;
                if(current_weight + upper > record_weight)
                {
                    expand(set2,i+1,upper);
                }
                free(set2);
                c[v] = record_weight;
            }
            j=0;
        }
nobs:
        /*
           Find the exact solution of the orignial input graph.
         */
        for(;i<number_of_subsets; i++)
        {
            for(;j<subset_size[i];j++)
            {
                set[i] += (1<<j);
            }
            j=0;
        }
        int upper=0;
        for(int i = 0; i < number_of_subsets; ++i)
        {
            upper += optimal_table[i][set[i]];
        }
        current_size = 0;
        current_weight = 0;
        if(current_weight + upper > record_weight)
        {
            expand(set,number_of_subsets,upper);
        }
    }

    free(set);
    free(c);
    free(current);
}

/*
   Branching procedure.
   <args>
    set: a vertex subset
    set_size: the size of "set"
    upper: an upper bound of the graph induced by "set"
 */
void expand(int *set,int set_size,int upper)
{
    ++branch_count;

    int i=set_size;
    while(i--)
    {
        while(set[i] != 0) 
        {
            /* check upper bound of optimal tables */
            if(current_weight + upper <= record_weight) 
            {
                return;
            }
            int msb=msb_table[set[i]];
            int vertex = (i*limit) + msb;
            /* check upper bound of c[]*/
            if(current_weight + c[vertex] <= record_weight)
            {
                return;
            }
            /* add vertex to current */
            current[current_size++] = vertex;
            current_weight += weight[vertex];
            /* make new set */
            int set2_size = (vertex-1)/limit+1;
            if(vertex==0)
            {
                set2_size=0;
            }
            int* set2 = (int *)malloc(sizeof(int)*set2_size);
            int *adjv = adjacency_matrix[vertex];
            int new_upper=0;
            {
                int j=set2_size;
                while(j--)
                {
                    set2[j] = set[j] & adjv[j];
                    new_upper+=optimal_table[j][set2[j]];
                }
            }
            if(current_weight + new_upper > record_weight)
            {
                expand(set2,set2_size,new_upper);
            }
            free(set2);
            --current_size;
            current_weight -= weight[vertex];
            /* delete vertex from set */
            upper-=optimal_table[i][set[i]];
            set[i] -= 1<<msb;
            upper+=optimal_table[i][set[i]];
        }
    }
    if(current_weight > record_weight) 
    {
        memcpy(record,current,sizeof(int)*current_size);
        record_size=current_size;
        record_weight = current_weight;
    }
}
