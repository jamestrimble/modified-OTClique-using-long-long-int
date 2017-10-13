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
#include "vertex_sequence_and_partition.h"
#include <stdlib.h>

static int* wt; /* vertex weight used in sorting*/
static int* dg; /* vertex degree used in sorting*/

/*
   comparing functions used in qsort()
 */
static int comp_weight_nondecreasing( const void *c1, const void *c2 )
{
    int v1=*(int *)c1;
    int v2=*(int *)c2;
    if(wt[v1] == wt[v2])
    {
        return dg[v2]-dg[v1]; /* degree nonincreasing */
    }
    return wt[v1]-wt[v2]; /* weight nondecreasing */
}

static int comp_degree_nondecreasing( const void *c1, const void *c2 )
{
    int v1=*(int *)c1;
    int v2=*(int *)c2;
    return dg[v1]-dg[v2]; /* degree nondecreasing */
}

/*
   make a vertex sequence and partition by greedy coloring for weighted case
 */
sequence_and_partition * coloring_weighted(weighted_graph *graph,int limit)
{
    int n=graph->n;
    int **adjacency_matrix=graph->adjacency_matrix;

    int color_size_limit=limit;
    {
        double edge_density=((double)graph->m/(n * (n-1) /2));
        if(edge_density>0.5) color_size_limit=limit;
        else if(edge_density>=0.4) color_size_limit=8;
        else if(edge_density>=0.3) color_size_limit=12;
        else if(edge_density>=0.2) color_size_limit=20;
        if(color_size_limit>limit) color_size_limit=limit;
    }

    /* sort vertices */
    wt=graph->weight;
    dg=(int *)calloc(n,sizeof(int));
    for(int i=0; i < n;i++) //calculate degree
    {
        int *adjv=adjacency_matrix[i];
        for(int j=0; j < n;j++)
        {
            if(adjv[j])
            {
                dg[i]++;
            }
        }
    }

    int *order=(int *)malloc(sizeof(int)*n);
    for(int i=0; i < n;i++)
    {
        order[i] = i;
    }
    qsort(order,n,sizeof(int),comp_weight_nondecreasing);

    /* create uncolored set (bit set) */
    int *uncolored=(int *)malloc(sizeof(int)*(n));
    for(int i=n-1;i>=0;i--)
    {
        uncolored[i]=1;
    }


    sequence_and_partition * result=(sequence_and_partition *)malloc(sizeof(sequence_and_partition));

    result->sequence=(int *)malloc(sizeof(int)*n);
    int *seq=result->sequence;

    /* greedy coloring */
    int number_of_colors=0;
    int *color_size=(int *)calloc(n,sizeof(int));
    {
        int i;
        for(int k = n; k>0;)
        {
            number_of_colors++;
            i=k;
            for(int j=n-1;j>=0;j--)
            {
                if(uncolored[j])
                {
                    int v=order[j];
                    int *adjv=adjacency_matrix[v];
                    int independent=1;
                    for(int h=i;h<k;++h)
                    {
                        if(adjv[seq[h]])
                        {
                            independent=0;
                            break;
                        }
                    }
                    if(independent)
                    {
                        seq[--i] = v;
                        uncolored[j]=0;
                        if(++color_size[number_of_colors-1] == color_size_limit)
                        {
                            goto OUT;
                        }
                    }
                }
            }
OUT:
            k=i;
        }
    }

    /* partition by colors */
    result->subset_size=(int *)calloc(number_of_colors,sizeof(int));
    int *size=result->subset_size;
    int num_of_subsets=1;
    for(int i=number_of_colors-1;i>=0;--i)
    {
        if(size[num_of_subsets-1]+color_size[i] > limit)
        {
            ++num_of_subsets;
        }
        size[num_of_subsets-1]+=color_size[i];
    }
    result->number_of_subsets=num_of_subsets;
    free(color_size);
    free(order);
    free(dg);
    free(uncolored);
    return result;
}

/*
   make a vertex sequence and partition by greedy coloring for unweighted case
 */
sequence_and_partition * coloring_unweighted(weighted_graph *graph,int limit)
{
    int n=graph->n;
    int **adjacency_matrix=graph->adjacency_matrix;

    int color_size_limit=limit;

    /* sort vertices */
    dg=(int *)calloc(n,sizeof(int));
    for(int i=0; i < n;i++) //calculate degree
    {
        int *adjv=adjacency_matrix[i];
        for(int j=0; j < n;j++)
        {
            if(adjv[j])
            {
                dg[i]++;
            }
        }
    }

    int *order=(int *)malloc(sizeof(int)*n);
    for(int i=0; i < n;i++)
    {
        order[i] = i;
    }
    qsort(order,n,sizeof(int),comp_degree_nondecreasing);

    /* create uncolored set (bit set) */
    int *uncolored=(int *)malloc(sizeof(int)*n);
    for(int i=n-1;i>=0;i--)
    {
        uncolored[i]=1;
    }


    sequence_and_partition * result=(sequence_and_partition *)malloc(sizeof(sequence_and_partition));

    result->sequence=(int *)malloc(sizeof(int)*n);
    int *seq=result->sequence;

    /* greedy coloring */

    int *color_size=(int *)calloc(n,sizeof(int));
    int number_of_colors=0;
    {
        int i;
        for(int k=n;k>0;)
        {
            number_of_colors++;
            i=k;
            for(int j=n-1;j>=0;j--)
            {
                if(uncolored[j])
                {
                    int v=order[j];
                    int *adjv=adjacency_matrix[v];
                    int independent=1;
                    for(int h=i;h<k;++h)
                    {
                        if(adjv[seq[h]])
                        {
                            independent=0;
                            break;
                        }
                    }
                    if(independent)
                    {
                        seq[--i] = v;
                        uncolored[j]=0;
                        if(++color_size[number_of_colors-1] ==color_size_limit)
                        {
                            goto OUT;
                        }
                    }
                }
            }
OUT:
            k=i;
        }
    }
    /* partition by colors */
    result->subset_size=(int *)calloc(number_of_colors,sizeof(int));
    int *size=result->subset_size;
    int num_of_subsets=1;
    for(int i=number_of_colors-1;i>=0;--i)
    {
        if(size[num_of_subsets-1]+color_size[i] > limit)
        {
            ++num_of_subsets;
        }
        size[num_of_subsets-1]+=color_size[i];
    }
    result->number_of_subsets=num_of_subsets;
    /*
       Reverse vertex sequence.
     */
    for(int i=0;i<num_of_subsets/2;i++)
    {
        int j=size[i];
        size[i]=size[num_of_subsets-1-i];
        size[num_of_subsets-1-i]=j;
    }
    for(int i=0;i<n/2;i++)
    {
        int j=seq[i];
        seq[i]=seq[n-1-i];
        seq[n-1-i]=j;
    }
    free(order);
    free(color_size);
    free(dg);
    free(uncolored);
    return result;
}
