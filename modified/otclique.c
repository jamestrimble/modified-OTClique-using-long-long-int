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

#include "c_program_timing.h"

static clock_t start_precomputation;
static clock_t end_precomuputation;
static clock_t start_branch_and_bound;
static clock_t end_branch_and_bound;
static unsigned long branch_count;
static long long **adjacency_matrix;
static long long *weight;
static long long *c;
static long long *record;
static long long record_weight;
static long long record_size;
static long long *current;
static long long current_size;
static long long current_weight;
static long long **optimal_table;
static long long limit;
static long long *msb_table;
static long long weighted;
static long long number_of_subsets,*subset_size;
static weighted_graph *input_graph;
static weighted_graph *reconstructed_graph;
static long long *seq=NULL;

static void precomputation();
static void branch_and_bound();
static void expand(long long *set,long long set_size,long long upper);

clique * otclique(weighted_graph *graph,long long subset_size_limit)
{
    double sec_precomputation;
    double sec_branch_and_bound;
    double sec_total;
    clique *maximum_weight_clique;
    long long n=graph->n;
    branch_count=0;
    input_graph=graph;
    limit=subset_size_limit;

    printf("Subset size limit = %lld\n", limit);
    start_precomputation=clock();

    /* check if the graph is weighted or unweighted */
    {
        weighted=0;
        long long weight0=graph->weight[0];
        long long *weight=graph->weight;
        for(long long i=0;i<n;i++)
        {
            if(weight0 != weight[i])
            {
                weighted=1;
                break;
            }
        }
    }

    /* create msb_table to get msb of any bit vector */
    msb_table = (long long *)malloc((1 << limit)*sizeof(long long));
    msb_table[0] = -1;
    for(long long i = 0; i < limit; ++i) 
    {
        long long from = 1 << i;
        long long to = 1 << (i + 1);
        for(long long j = from; j < to; ++j) 
        {
            msb_table[j] = i;
        }
    }

    /* precomputation phase */
    precomputation();

    end_precomuputation=clock();
    sec_precomputation=(double)(end_precomuputation-start_precomputation)/CLOCKS_PER_SEC;
    printf("%lld subsets created from %lld vertices \n", number_of_subsets, n);
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
    maximum_weight_clique->set=(long long *)malloc(sizeof(long long) * record_size);
    for(long long i = 0; i < record_size; ++i)
    {
        maximum_weight_clique->set[i]=seq[record[i]];
    }

    for(long long i = 0;i<number_of_subsets;i++)
    {
        free(optimal_table[i]);
    }
    free(optimal_table);
    free(msb_table);
    {
        long long n_r=reconstructed_graph->n;
        for(long long i=0; i<n_r; i++)
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

    long long *seq0=seq_and_partition->sequence;
    number_of_subsets=seq_and_partition->number_of_subsets;
    subset_size = seq_and_partition->subset_size;
    free(seq_and_partition);

    /* encode to bit vector*/
    seq=(long long *)calloc(number_of_subsets*limit,sizeof(long long));
    {
        long long k=0;
        for(long long i=0;i<number_of_subsets;i++)
        {
            long long h=i*limit;
            for(long long l=0;l<subset_size[i];l++)
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
    long long n=input_graph->n;

    /* initialize variables used in branch-and-bound phase */
    record=(long long *)malloc(sizeof(long long)*n);
    record_weight=0;
    record_size=0;
    current = (long long *)malloc(sizeof(long long)*n);
    c = (long long *)malloc(sizeof(long long) * (number_of_subsets*limit));
    for(long long i=0;i<number_of_subsets*limit;i++)
    {
        c[i]=LLONG_MAX/2;
    }

    long long stop=n;
    if(weighted)
    {
        stop=n*0.8;
    }

    /* main loop */
    long long *set = (long long *)calloc(number_of_subsets,sizeof(long long));
    {
        long long i=0;long long j=0;long long l=0;
        for(i=0; i<number_of_subsets; i++)
        {
            for(;j<subset_size[i];j++)
            {
                if(l++==stop)
                {
                    goto nobs; /* stop calculation of c[] */
                }
                set[i] += (1<<j);
                long long v = (i*limit) + j;
                long long* adjv=adjacency_matrix[v];
                long long* set2=(long long *)calloc(i+1,sizeof(long long));
                long long upper=0;
                long long k=(v-1)/limit+1;
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
        long long upper=0;
        for(long long i = 0; i < number_of_subsets; ++i)
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
void expand(long long *set,long long set_size,long long upper)
{
    ++branch_count;
    if (branch_count % 100000 == 0)
        check_for_timeout();
    if (is_timeout_flag_set())
        return;

    long long i=set_size;
    while(i--)
    {
        while(set[i] != 0) 
        {
            /* check upper bound of optimal tables */
            if(current_weight + upper <= record_weight) 
            {
                return;
            }
            long long msb=msb_table[set[i]];
            long long vertex = (i*limit) + msb;
            /* check upper bound of c[]*/
            if(current_weight + c[vertex] <= record_weight)
            {
                return;
            }
            /* add vertex to current */
            current[current_size++] = vertex;
            current_weight += weight[vertex];
            /* make new set */
            long long set2_size = (vertex-1)/limit+1;
            if(vertex==0)
            {
                set2_size=0;
            }
            long long* set2 = (long long *)malloc(sizeof(long long)*set2_size);
            long long *adjv = adjacency_matrix[vertex];
            long long new_upper=0;
            {
                long long j=set2_size;
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
        memcpy(record,current,sizeof(long long)*current_size);
        record_size=current_size;
        record_weight = current_weight;
    }
}
