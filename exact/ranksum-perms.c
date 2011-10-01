#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

static inline int rand_int(int n) {
  int limit = RAND_MAX - RAND_MAX % n;
  int rnd;
  do {
    rnd = rand();
  }
  while (rnd >= limit);
  return rnd % n;
}

void shuffle(int *annot, int n) {
  for (int i = n - 1; i > 0; i--) {
    int j = rand_int(i + 1);
    int temp = annot[i];
    annot[i] = annot[j];
    annot[j] = temp;
  }
}

static inline void update(int *annot, int *rank, long *result) {
  if (*annot == 1) {
    *result += *rank;
  }
}

long ranksum(int *annots, int n, int offset) {
  long result = 0;
  int rank = 0;
  for (int *ap = annots + offset; ap < annots + n; ap++, rank++) {
    update(ap, &rank, &result);
  }
  for (int *ap = annots; ap < annots + offset; ap++, rank++) {
    update(ap, &rank, &result);
  }
  fprintf(stderr, "%ld\n", result);
  return result;
}

int main(int argc, char **argv) {
  srand(time(NULL));
  int n;
  scanf("%d", &n);
  int *annots = calloc(n, sizeof(int));
  assert(annots != NULL);
  for(int *ap = annots; ap < annots + n; ap++) {
    scanf("%d", ap);
  }
  long orig = ranksum(annots, n, 0);
  int count = 0;
  for (int i = 1; i < n; i++) {
    shuffle(annots, n);
    if (ranksum(annots, n, 0) <= orig) {
      count++;
    }
  }
  printf("%.3e\n", (double)count / n);
}
