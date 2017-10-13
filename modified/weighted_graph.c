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
#include <stdio.h>
#include <stdlib.h>
#include <memory.h>

/*
   Read a weighted DIMACS format graph.
   <args>
    inFile: input graph file name
   <return>
    a graph or NULL(can't read file)
 */
weighted_graph * read_graph(char *inFile)
{
    weighted_graph *graph;
    FILE *fp=fopen(inFile,"r");
    if(fp == NULL)
    {
        fprintf(stderr,"Can't read \"%s\"\n",inFile);
        exit(1);
    }
    else
    {
        graph=(weighted_graph *)malloc(sizeof(weighted_graph));
        long long n;
        long long **adjacency_matrix=NULL;
        long long *weight=NULL;
        char form;
        while(fscanf(fp," %c",&form) != EOF)
        {
            switch(form)
            {
                case 'c':
                    if( fscanf(fp,"%*[^\n]") != 0 )
                    {
                        fprintf(stderr,"input file error\n");
                        exit(1);
                    }
                    break;
                case 'd':
                    if( fscanf(fp,"%*[^\n]") != 0 )
                    {
                        fprintf(stderr,"input file error\n");
                        exit(1);
                    }
                    break;
                case 'v':
                    if( fscanf(fp,"%*[^\n]") != 0 )
                    {
                        fprintf(stderr,"input file error\n");
                        exit(1);
                    }
                    break;
                case 'x':
                    if( fscanf(fp,"%*[^\n]") != 0 )
                    {
                        fprintf(stderr,"input file error\n");
                        exit(1);
                    }
                    break;
                case 'p':
                    if( fscanf(fp,"%*s %lld %lld",&graph->n,&graph->m) != 2 )
                    {
                        fprintf(stderr,"input file error\n");
                        exit(1);
                    }
                    n=graph->n;
                    graph->weight=(long long *)malloc(n*sizeof(long long));
                    weight=graph->weight;
                    for(long long i=0;i<n;i++)
                    {
                        weight[i]=1;
                    }
                    graph->adjacency_matrix=(long long **)malloc(n * sizeof(long long *) + n * n * sizeof(long long));
                    adjacency_matrix=graph->adjacency_matrix;
                    adjacency_matrix[0]=(long long *)(adjacency_matrix + n);
                    for(long long i=1;i<n;i++)
                    {
                        adjacency_matrix[i] = adjacency_matrix[0] + i * n;
                    }
                    for(long long i=0;i<n;i++)
                    {
                        for(long long j=0;j<n;j++)
                        {
                            adjacency_matrix[i][j] = 0;
                        }
                    }
                    break;
                case 'e':
                    {
                        long long e1,e2;
                        if( fscanf(fp,"%lld %lld",&e1,&e2) != 2 )
                        {
                            fprintf(stderr,"input file error\n");
                            exit(1);
                        }
                        adjacency_matrix[e1-1][e2-1]=1;
                        adjacency_matrix[e2-1][e1-1]=1;
                    }
                    break;
                case 'n':
                    {
                        long long v,w;
                        if( fscanf(fp,"%lld %lld",&v,&w) != 2 )
                        {
                            fprintf(stderr,"input file error\n");
                            exit(1);
                        }
                        weight[v-1]=w;
                    }
                    break;
                default:
                    fprintf(stderr,"input file error\n");
                    exit(1);
                    break;
            }
        }
    }
    fclose(fp);
    return graph;
}

/*
   Create a subgraph induced by the given vertex set.
   <args>
    seq: vertex subset of graph
    size: number of vertices in seq
    graph: a graph
   <return>
    vertex induced subgraph
 */
weighted_graph * create_vertex_induced_subgraph(long long *seq,long long size,weighted_graph *graph)
{
    long long *graph_weight=graph->weight;
    long long **graph_adjacency_matrix=graph->adjacency_matrix;

    weighted_graph *vertex_induced_subgraph=(weighted_graph *)malloc(sizeof(weighted_graph));

    /* initialize */
    vertex_induced_subgraph->n=size;
    vertex_induced_subgraph->adjacency_matrix=(long long **)malloc(size * sizeof(long long *) + size * size * sizeof(long long));
    vertex_induced_subgraph->weight = (long long *)malloc(sizeof(long long) * vertex_induced_subgraph->n);

    /* create adjacency matrix */
    {
        long long **subgraph_adjacency_matrix=vertex_induced_subgraph->adjacency_matrix;
        subgraph_adjacency_matrix[0]=(long long *)(subgraph_adjacency_matrix + size);
        for(long long i=0;i<size;i++)
        {
            subgraph_adjacency_matrix[i] = subgraph_adjacency_matrix[0] + i * size;

            long long* adj=graph_adjacency_matrix[seq[i]];
            for(long long j=0;j<size;j++)
            {
                subgraph_adjacency_matrix[i][j] = adj[seq[j]];
            }
        }
    }

    /* copy vertex weight */
    {
        long long *subgraph_weight=vertex_induced_subgraph->weight;
        for(long long i=0;i<size;i++)
        {
            subgraph_weight[i]=graph_weight[seq[i]];
        }
    }
    return vertex_induced_subgraph;
}

/*
   print the graph
   <args>
    graph: a graph
 */
void print_graph(weighted_graph *graph)
{
    long long n=graph->n;
    long long *weight=graph->weight;
    long long **adjacency_matrix=graph->adjacency_matrix;

    printf("%lld\n",n);

    for(long long i=0;i<n;i++)
    {
        printf("%lld ",weight[i]);
        long long *adj=adjacency_matrix[i];
        for(long long j=0;j<n;j++)
        {
            if(adj[j])
            {
                printf(" %lld",j);
            }
        }
        printf("\n");
    }
}

/*
   Get the adjacency matrix implemented by bit vector.
   (only bottom triangle)
   <args>
    graph: a graph
    unit: the length of one word bit vector
   <return>
    pointer to the adjacency matrix implemented by bit vector.
 */
long long ** get_bit_vector_adjacency_matrix(weighted_graph *graph,long long unit)
{
    long long n=graph->n;
    long long **adjacency_matrix=graph->adjacency_matrix;

    long long **bit_adj = (long long **)malloc(sizeof(long long *) * n);
    bit_adj[0]=(long long *)calloc(1,sizeof(long long));

    for(long long i = 1; i < n; i++) 
    {
        long long len = (i-1)/unit+1;
        bit_adj[i]= (long long *)calloc(len,sizeof(long long));
        long long *adji=adjacency_matrix[i];
        for(long long j=0;j<i;j++)
        {
            if(adji[j])
            {
                bit_adj[i][j/unit] += 1<<(j%unit);
            }
        }
    }
    return bit_adj;
}

/*
   get complement graph.
   <args>
    graph: input graph to make complement graph
   <return>
    complement graph
 */
weighted_graph * get_complement_graph(weighted_graph *graph)
{
    long long n=graph->n;

    weighted_graph *complement_graph=(weighted_graph *)malloc(sizeof(weighted_graph));

    complement_graph->n=n;
    complement_graph->m=(n * (n-1) /2)-graph->m;
    complement_graph->weight=(long long *)malloc(n*sizeof(long long));
    
    /* copy weight */
    memcpy(complement_graph->weight,graph->weight,sizeof(long long)*n);

    /* create adjacency_matrix */
    complement_graph->adjacency_matrix=(long long **)malloc(n * sizeof(long long *) + n * n * sizeof(long long));

    long long **adjacency_matrix=complement_graph->adjacency_matrix;
    adjacency_matrix[0]=(long long *)(adjacency_matrix + n);

    for(long long i=0;i<n;i++)
    {
        adjacency_matrix[i] = adjacency_matrix[0] + i * n;
        long long *adji=adjacency_matrix[i];
        long long *adji_origin=graph->adjacency_matrix[i];
        for(long long j=0;j<n;j++)
        {
            if( (!adji_origin[j]) && (i!=j) )
            {
                adji[j] = 1;
            }
            else
            {
                adji[j] = 0;
            }
        }
    }
    return complement_graph;
}
