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
        int n;
        int **adjacency_matrix=NULL;
        int *weight=NULL;
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
                    if( fscanf(fp,"%*s %d %d",&graph->n,&graph->m) != 2 )
                    {
                        fprintf(stderr,"input file error\n");
                        exit(1);
                    }
                    n=graph->n;
                    graph->weight=(int *)malloc(n*sizeof(int));
                    weight=graph->weight;
                    for(int i=0;i<n;i++)
                    {
                        weight[i]=1;
                    }
                    graph->adjacency_matrix=(int **)malloc(n * sizeof(int *) + n * n * sizeof(int));
                    adjacency_matrix=graph->adjacency_matrix;
                    adjacency_matrix[0]=(int *)(adjacency_matrix + n);
                    for(int i=1;i<n;i++)
                    {
                        adjacency_matrix[i] = adjacency_matrix[0] + i * n;
                    }
                    for(int i=0;i<n;i++)
                    {
                        for(int j=0;j<n;j++)
                        {
                            adjacency_matrix[i][j] = 0;
                        }
                    }
                    break;
                case 'e':
                    {
                        int e1,e2;
                        if( fscanf(fp,"%d %d",&e1,&e2) != 2 )
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
                        int v,w;
                        if( fscanf(fp,"%d %d",&v,&w) != 2 )
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
weighted_graph * create_vertex_induced_subgraph(int *seq,int size,weighted_graph *graph)
{
    int *graph_weight=graph->weight;
    int **graph_adjacency_matrix=graph->adjacency_matrix;

    weighted_graph *vertex_induced_subgraph=(weighted_graph *)malloc(sizeof(weighted_graph));

    /* initialize */
    vertex_induced_subgraph->n=size;
    vertex_induced_subgraph->adjacency_matrix=(int **)malloc(size * sizeof(int *) + size * size * sizeof(int));
    vertex_induced_subgraph->weight = (int *)malloc(sizeof(int) * vertex_induced_subgraph->n);

    /* create adjacency matrix */
    {
        int **subgraph_adjacency_matrix=vertex_induced_subgraph->adjacency_matrix;
        subgraph_adjacency_matrix[0]=(int *)(subgraph_adjacency_matrix + size);
        for(int i=0;i<size;i++)
        {
            subgraph_adjacency_matrix[i] = subgraph_adjacency_matrix[0] + i * size;

            int* adj=graph_adjacency_matrix[seq[i]];
            for(int j=0;j<size;j++)
            {
                subgraph_adjacency_matrix[i][j] = adj[seq[j]];
            }
        }
    }

    /* copy vertex weight */
    {
        int *subgraph_weight=vertex_induced_subgraph->weight;
        for(int i=0;i<size;i++)
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
    int n=graph->n;
    int *weight=graph->weight;
    int **adjacency_matrix=graph->adjacency_matrix;

    printf("%d\n",n);

    for(int i=0;i<n;i++)
    {
        printf("%d ",weight[i]);
        int *adj=adjacency_matrix[i];
        for(int j=0;j<n;j++)
        {
            if(adj[j])
            {
                printf(" %d",j);
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
int ** get_bit_vector_adjacency_matrix(weighted_graph *graph,int unit)
{
    int n=graph->n;
    int **adjacency_matrix=graph->adjacency_matrix;

    int **bit_adj = (int **)malloc(sizeof(int *) * n);
    bit_adj[0]=(int *)calloc(1,sizeof(int));

    for(int i = 1; i < n; i++) 
    {
        int len = (i-1)/unit+1;
        bit_adj[i]= (int *)calloc(len,sizeof(int));
        int *adji=adjacency_matrix[i];
        for(int j=0;j<i;j++)
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
    int n=graph->n;

    weighted_graph *complement_graph=(weighted_graph *)malloc(sizeof(weighted_graph));

    complement_graph->n=n;
    complement_graph->m=(n * (n-1) /2)-graph->m;
    complement_graph->weight=(int *)malloc(n*sizeof(int));
    
    /* copy weight */
    memcpy(complement_graph->weight,graph->weight,sizeof(int)*n);

    /* create adjacency_matrix */
    complement_graph->adjacency_matrix=(int **)malloc(n * sizeof(int *) + n * n * sizeof(int));

    int **adjacency_matrix=complement_graph->adjacency_matrix;
    adjacency_matrix[0]=(int *)(adjacency_matrix + n);

    for(int i=0;i<n;i++)
    {
        adjacency_matrix[i] = adjacency_matrix[0] + i * n;
        int *adji=adjacency_matrix[i];
        int *adji_origin=graph->adjacency_matrix[i];
        for(int j=0;j<n;j++)
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
