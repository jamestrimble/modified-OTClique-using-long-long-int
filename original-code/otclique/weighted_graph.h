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

#ifndef weighted_graph_h
#define weighted_graph_h

typedef struct
{
  int n; /* number of vertices */
  int m; /* number of edges */
  int **adjacency_matrix;
  int *weight;
} weighted_graph;

weighted_graph * read_graph(char *inFile);
weighted_graph * create_vertex_induced_subgraph(int *seq,int size,weighted_graph *graph);
void print_graph(weighted_graph *graph);
int ** get_bit_vector_adjacency_matrix(weighted_graph *graph,int unit);
weighted_graph * get_complement_graph(weighted_graph *graph);

#endif
