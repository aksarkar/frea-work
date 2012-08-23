#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int cmp(const void *l, const void *r) {
  double x = *(double*)l;
  double y = *(double*)r; 
  if (x < y) {
    return -1;
  }
  else if (x > y) {
    return 1;
  }
  else {
    return 0;
  }
}

static inline void update(int *annot, int *rank, long *result) {
  if (*annot == 1) {
    *result += *rank;
  }
}

long ranksum(int *rs, int *xs, int n, int offset) {
  long result = 0;
  int *r = rs;
  for (int *ap = xs + offset; ap < xs + n; ap++, r++) {
    update(ap, r, &result);
  }
  for (int *ap = xs; ap < xs + offset; ap++, r++) {
    update(ap, r, &result);
  }
  return result;
}

int main(int argc, char **argv) {
  int n;
  scanf("%d", &n);
  double *ps = calloc(n, sizeof(double));
  assert(ps != NULL);
  int *xs = calloc(n, sizeof(int));
  assert(xs != NULL);
  {
    double *p;
    int *x;
    for(p = ps, x = xs; x < xs + n; p++, x++) {
      scanf("%lf %d", p, x);
    }
  }

  double *qs = calloc(n, sizeof(double));
  assert(qs != NULL);
  memcpy(qs, ps, n * sizeof(double));
  qsort(qs, n, sizeof(double), cmp);
  int *rs = calloc(n, sizeof(int));
  assert(rs != NULL);
  {
    double *p;
    int *r;
    for (p = ps, r = rs; p < ps + n; p++, r++) {
      void *t = bsearch(p, qs, n, sizeof(double), cmp);
      assert(t != NULL);
      *r = (double*)t - qs;
    }
  }
  free(ps);
  free(qs);

  long orig = ranksum(rs, xs, n, 0);
  int count = 0;
  for (int i = 1; i < n; i++) {
    if (ranksum(rs, xs, n, i) <= orig) {
      count++;
    }
  }
  printf("%.3e\n", (double)count / n);
}
