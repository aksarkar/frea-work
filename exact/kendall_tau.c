#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

/**
   AVL tree
 */
typedef struct avl {
  double key;
  struct avl *left, *right;
  int height;

  // needed for Kendall tau
  int count;
  int nbefore;
} avl_t;

static inline int max(int a, int b) {
  return a > b ? a : b;
}

static inline int avl_count(avl_t *avl) {
  return avl == NULL ? 0 : avl->count;
}

static inline int avl_height(avl_t *avl) {
  return avl == NULL ? -1 : avl->height;
}

static inline void avl_update(avl_t *avl) {
  avl->count = avl_count(avl->right) + avl_count(avl->left) + 1;
  avl->height = max(avl_height(avl->left), avl_height(avl->right)) + 1;
}

avl_t *avl_rotl(avl_t *avl) {
  assert(avl != NULL);
  assert(avl->left != NULL);
  avl_t *result = avl->left;
  avl->left = result->right;
  result->right = avl;
  avl_update(avl);
  avl_update(result);
  return result;
}

avl_t *avl_rotr(avl_t *avl) {
  assert(avl != NULL);
  assert(avl->right != NULL);
  avl_t *result = avl->right;
  avl->right = result->left;
  result->left = avl;
  avl_update(avl);
  avl_update(result);
  return result;
}

avl_t *avl_insert(avl_t *avl, double key, int *nbefore) {
  if (avl == NULL) {
    avl = malloc(sizeof(avl_t));
    assert(avl != NULL);
    avl->key = key;
    avl->left = NULL;
    avl->right = NULL;
    avl->height = 0;
    avl->count = 0;
  }
  else if (key < avl->key) {
    avl->left = avl_insert(avl->left, key, nbefore);
    if (avl_height(avl->left) > avl_height(avl->right) + 1) {
      if (key < avl->left->key) {
        avl = avl_rotl(avl);
      }
      else {
        avl->left = avl_rotr(avl->left);
        avl = avl_rotl(avl);
      }
    }
  }
  else if (key > avl->key) {
    *nbefore += avl_count(avl->left) + 1;
    avl->right = avl_insert(avl->right, key, nbefore);
    if (avl_height(avl->right) > avl_height(avl->left) + 1) {
      if (key > avl->right->key) {
        avl = avl_rotr(avl);
      }
      else {
        avl->right = avl_rotl(avl->right);
        avl = avl_rotr(avl);
      }
    }
  }
  avl_update(avl);
  return avl;
}

void avl_free(avl_t *avl) {
  if (avl == NULL) {
    return;
  }
  avl->left != NULL ? avl_free(avl->left) : 0;
  avl->right != NULL ? avl_free(avl->right) : 0;
  free(avl);
}

void nbefore(double *ys, int *zs, int n) {
  avl_t *avl = NULL;
  for (int i = 0; i < n; i++) {
    avl = avl_insert(avl, ys[i], &zs[i]);
  }
  free(avl);
}

static inline void update(double *ys, int *zs, int i, int j, long *test) {
  if (ys[i] > ys[j]) {
    zs[j] = max(zs[j] - 1, 0);
    *test += zs[j];
  }
  else {
    zs[i]++;
    *test++;
  }
}

int main() {
  int n;
  scanf("%d\n", &n);
  double *ys = calloc(n, sizeof(double));
  assert(ys != NULL);
  for (int i = 0; i < n; i++) {
    scanf("%lf", &ys[i]);
  }
  int *zs = calloc(n, sizeof(int));
  nbefore(ys, zs, n);
  long orig = 0;
  for (int i = 0; i < n; i++) {
    orig += zs[i];
  }
  int count = 0;
  for (int i = 1; i < n; i++) {
    long test = 0;
    zs[i] = 0;
    for (int j = i + 1; j < n; j++) {
      update(ys, zs, i, j, &test);
    }
    for (int j = 0; j < i; j++) {
      update(ys, zs, i, j, &test);
    }
    if (test >= orig) {
      count++;
    }
  }
  printf("%d\n", count);
}
