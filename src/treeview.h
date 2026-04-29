/*
 * treeview : Auxiliary routines for printing a IMG_NODE_t datastructure
 *
 * Adapted from the implementation documented originally on a blog at openasthra.com:
 *
 *    http://web.archive.org/web/20090617110918/http://www.openasthra.com/c-tidbits/printing-binary-trees-in-ascii/
 *
 * For now this is implemented completely as a self contained header file.
 *
 * -RMH
 */
#ifndef __TREEVIEW_H__
#define __TREEVIEW_H__


#include "msps.h"
#include "eleredef.h"


#define MAX_HEIGHT 1000
#define INFIN (1<<20)

typedef struct asciinode_struct asciinode;

struct asciinode_struct
{
  asciinode * left, * right;

  //length of the edge from this node to its children
  int edge_length;

  int height;

  int lablen;

  //-1=I am left, 0=I am root, 1=right
  int parent_dir;

  //max supported unit32 in dec, 10 digits max
  char label[20];
};


/* function prototypes -- implementations are in treeview.c */
void     print_ascii_tree(IMG_NODE_t *t);
void     print_level(int *print_next, asciinode *node, int x, int level);
void     compute_edge_lengths(int gap, int *lprofile, int *rprofile, asciinode *node);
asciinode *build_ascii_tree_recursive(IMG_NODE_t *t);
asciinode *build_ascii_tree(IMG_NODE_t *t);
void     free_ascii_tree(asciinode *node);
void     compute_lprofile(int *lprofile, asciinode *node, int x, int y);
void     compute_rprofile(int *rprofile, asciinode *node, int x, int y);
int      MIN(int X, int Y);
int      MAX(int X, int Y);

#endif /* __TREEVIEW_H__ */
