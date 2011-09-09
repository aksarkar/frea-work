#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

int ranksum(int *annots, int n, int offset) {
  int result = 0;
  int rank = 0;
  for (int *ap = annots + offset; ap < annots + n; ap++, rank++) {
    if (*ap == 1) {
      result += rank;
    }
  }
  for (int *ap = annots; ap < annots + offset; ap++, rank++) {
    if (*ap == 1) {
      result += rank;
    }
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
  int orig = ranksum(annots, n, 0);
  int count = 0;
  for (int i = 1; i < n; i++) {
    if (ranksum(annots, n, i) <= orig) {
      count++;
    }
  }
  printf("%.3e\n", (double)count / n);
}
