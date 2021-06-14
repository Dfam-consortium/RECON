typedef struct img_node {
  short recorded;
  IMAGE_t *to_image;
  struct img_node *sib;
  struct img_node *children;
} IMG_NODE_t;

