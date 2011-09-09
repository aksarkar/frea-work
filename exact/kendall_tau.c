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
  int count;  // needed for Kendall tau
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
  avl->count = avl_count(avl->left) + avl_count(avl->right) + 1;
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

avl_t *avl_insert(avl_t *avl, double key) {
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
    avl->left = avl_insert(avl->left, key);
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
    avl->right = avl_insert(avl->right, key);
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

int avl_lookup(avl_t *avl, double key) {
  if (avl == NULL) {
    return 0;
  }
  if (avl->key == key) {
    return avl_count(avl->left);
  }
  else if (avl->key > key) {
    return avl_lookup(avl->left, key);
  }
  else {
    return 1 + avl_count(avl->left) + avl_lookup(avl->right, key);
  }
}

void avl_free(avl_t *avl) {
  if (avl == NULL) {
    return;
  }
  avl->left != NULL ? avl_free(avl->left) : 0;
  avl->right != NULL ? avl_free(avl->right) : 0;
  free(avl);
}

double kendall_tau(double *ys, int n, int offset) {
  int count = 0;
  avl_t *avl = NULL;
  for (int i = offset; i < n; i++) {
    avl = avl_insert(avl, ys[i]);
    count += avl_lookup(avl, ys[i]);
  }
  for (int i = 0; i < offset; i++) {
    avl = avl_insert(avl, ys[i]);
    count += avl_lookup(avl, ys[i]);
  }
  free(avl);
  return 4. * count / (n * (n - 1)) - 1;
}

int main() {
  int n;
  scanf("%d\n", &n);
  double *ys = calloc(n, sizeof(double));
  assert(ys != NULL);
  for (int i = 0; i < n; i++) {
    scanf("%lf", &ys[i]);
  }
  printf("%lf\n", kendall_tau(ys, n, 0));
}
