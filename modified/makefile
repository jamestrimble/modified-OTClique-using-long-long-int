CC = gcc
CFLAGS = -O2 -Wall -pedantic-errors -std=c11
LDFLAGS =
.c.o:
	${CC} ${CFLAGS} -c $<

.PHONY : all
all: otclique otclique_mwvc

otclique : mwc.o otclique.o optimal_table.o vertex_sequence_and_partition.o clique.o weighted_graph.o
	${CC} ${CFLAGS} -o $@ $^ ${LDFLAGS}

otclique_mwvc : mwvc.o otclique.o optimal_table.o vertex_sequence_and_partition.o clique.o weighted_graph.o
	${CC} ${CFLAGS} -o $@ $^ ${LDFLAGS}

.PHONY : clean
clean:
	rm -f otclique otclique_mwvc *.o
