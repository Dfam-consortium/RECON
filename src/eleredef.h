#ifndef __ELEREDEF_H__
#define __ELEREDEF_H__

typedef struct img_node {
  short recorded;
  IMAGE_t *to_image;
  struct img_node *sib;
  struct img_node *children;
} IMG_NODE_t;

#endif
