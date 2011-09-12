#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

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
  return result;
}

int main(int argc, char **argv) {
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
    if (ranksum(annots, n, i) <= orig) {
      count++;
    }
  }
  printf("%.3e\n", (double)count / n);
}
