#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

int ranksum(char *annots, int n, int offset) {
  int result = 0;
  int rank = 0;
  for (char *ap = annots + offset; ap < annots + n; ap++, rank++) {
    if (*ap == '1') {
      result += rank;
    }
  }
  for (char *ap = annots; ap < annots + offset; ap++, rank++) {
    if (*ap == '1') {
      result += rank;
    }
  }
  return result;
}

int main(int argc, char **argv) {
  int n;
  scanf("%d", &n);
  char *annots = malloc(n);
  assert(annots != NULL);
  for(char *ap = annots; ap < annots + n; ap++) {
    scanf("%c", ap);
  }
  fclose(stdin);
  int orig = ranksum(annots, n, 0);
  int count = 0;
  for (int i = 1; i < n; i++) {
    if (ranksum(annots, n, i) <= orig) {
      count++;
    }
  }
  printf("%.3e\n", (double)count / n);
}
