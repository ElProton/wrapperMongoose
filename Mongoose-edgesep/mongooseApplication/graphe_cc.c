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
	
	int k = 1;  /* number of connected components in the complement graph */
	int bot = 0;
	int mid = 0;
	
	/* Visited vertices are in queue[0:mid]. 
	   The current BFS queue is queue[bot:mid]
	   Unvisited vertices are in queue[mid:n]. 
	*/

	for (int i = 0; i < n; i++) {
		assert(mid == bot);

		if (bot == n)
			break;

		mid++;

		while (bot < mid && mid <= n) {
			int u = queue[bot++];
			mark[u] = k;

			for (int i = 0; i < n; i++)
				printf("%d : %d\n",i,mark[i]);
			
			for (int it = Ap[u]; it < Ap[u + 1]; it++){
				
				if (mark[Aj[it]] == 0) {
					mark[Aj[it]] = k;
					spasm_swap(queue, mid, Aj[it]);
					mid++;
					printf("swap %d to %d\n",mid,Aj[it]);
				}
			}
			printf("bot : %d, mid : %d\n",bot,mid);
		}

		k++;

		bot = mid;
	}
	
	/* just print the number of components */
	printf("%d\n", k);
	
	free(queue);
	free(mark);
	spasm_csr_free(A);
	return 0;
}
