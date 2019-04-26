#include <assert.h>
#include <stdio.h>

#include "mini_spasm.h"

/* compute the connected components of the complement graph. */


int main() 
{
	spasm_triplet * T = spasm_load_mm(stdin, -1);
	assert (T->n == T->m);

	spasm * A = spasm_compress(T);
	spasm_triplet_free(T);
	int n = A->n;
	int *Ap = A->p;
	int *Aj = A->j;

	/**************** BFS ****************/
	int *queue = spasm_malloc(sizeof(int) * n);
	int *mark = spasm_malloc(sizeof(int) * n);


	/* there are two data structures:
	   - a set of unvisited vertices
	   - a BFS queue of visited vertices
	*/
	for (int i = 0; i < n; i++)
		queue[i] = i;
	for (int i = 0; i < n; i++)
		mark[i] = 0;
	
	int k = 0;  /* number of connected components in the complement graph */
	int bot = 0;
	int mid = 0;
	
	/* Visited vertices are in queue[0:mid]. 
	   The current BFS queue is queue[bot:mid]
	   Unvisited vertices are in queue[mid:n]. 
	*/

	for (int i = 0; i < n; i++) {
		assert(mid == bot);

		if (bot == n)
			break;  /* all vertices have been visited. */
		
		/* start BFS in the complement graph starting at the next unvisited vertex. 
		   Move it from the unvisited list to the current BFS queue. 
		   We compute its connected component in the complement graph. */
		mid++;

		/* while the BFS queue is not empty */
		while (bot < mid) {
			int u = queue[bot++];  /* dequeue a vertex  */
			
			/* mark all its neighbors -- they are non-neighbors in the complement */
			for (int it = Ap[u]; it < Ap[u + 1]; it++)
				mark[Aj[it]] = 1;
			
			/* examine all unvisited vertices. If they are marked (neighbors of u in G),
			   they are not adjacent to u in the complement graph, so they remain unvisited.
			   Otherwise, we move them to the current BFS queue because they belong to the 
			   connected component in the complement graph and we had not seen them before.
			   When this is over, queue[hi:n] are unvisited, while queue[mid:hi] have been 
			   visited from u */
			int hi = n;
			while (mid < hi)
				if (mark[queue[mid]]) {    /* A: non-neighbor: remain unvisited */
					hi--;
					spasm_swap(queue, mid, hi);
				} else {
					mid++;             /* B: neighbor: move to current BFS queue */
				}

			/* erase all marks. */
			for (int it = Ap[u]; it < Ap[u + 1]; it++)
				mark[Aj[it]] = 0;

			/* complexity analysis: 
			     A can happen at most once per edge during the whole algorithm.
			     B can happen at most once per vertex (once visited, stay visited).
			     updates to the mark array happen once per edge during the whole algorithm.
			     ---> the whole thing is linear
			*/

			for (int i = 0; i < n; i++)
				printf("%d : %d\n",i,mark[i]);
		}

		/* BFS is over. This makes one more connected component in the complement graph.
		   Its vertices are in queue[bot:mid]. */
		k++;

		/* ready for the next one. Empty the BFS queue */
		bot = mid;
	}
	
	/* just print the number of components */
	printf("%d\n", k);
	
	free(queue);
	free(mark);
	spasm_csr_free(A);
	return 0;
}
