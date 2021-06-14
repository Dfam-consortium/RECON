/*
The below routine prints tree in ascii for a given Tree representation which contains list of nodes, and node structure is this

    struct Tree
    {
      Tree * left, * right;
      int element;
    };

This pic illustrates what the below routine does on canvas..
ascii tree

Here is the printing routine..

  http://web.archive.org/web/20090617110918/http://www.openasthra.com/c-tidbits/printing-binary-trees-in-ascii/

Auxiliary routines..
*/

//#include <stdio.h>
//#include <stdlib.h>
//#include <string.h>

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


void print_ascii_tree(IMG_NODE_t * t);
void print_level(int *printnext, asciinode *node, int x, int level);
void compute_edge_lengths(int gap, int *lprofile, int *rprofile, asciinode *node);
asciinode * build_ascii_tree_recursive(IMG_NODE_t * t);
asciinode * build_ascii_tree(IMG_NODE_t * t);
void free_ascii_tree(asciinode *node);
void compute_lprofile(int *lprofile, asciinode *node, int x, int y);
void compute_rprofile(int *rprofile, asciinode *node, int x, int y);
int MIN (int X, int Y);
int MAX (int X, int Y);


////////////////////////////////////

//prints ascii tree for given Tree structure
void print_ascii_tree(IMG_NODE_t * t)
{
  int lprofile[MAX_HEIGHT];
  int rprofile[MAX_HEIGHT];
  //used for printing next node in the same level,
  //this is the x coordinate of the next char printed
  int print_next = 0;
  //adjust gap between left and right nodes
  int gap = 7;
  asciinode *proot;
  int xmin, i;
  if (t == NULL) return;
  proot = build_ascii_tree(t);
  compute_edge_lengths(gap, lprofile, rprofile, proot);
  for (i=0; i<proot->height && i < MAX_HEIGHT; i++)
  {
    lprofile[i] = INFIN;
  }
  compute_lprofile(lprofile, proot, 0, 0);
  xmin = 0;
  for (i = 0; i < proot->height && i < MAX_HEIGHT; i++)
  {
    xmin = MIN(xmin, lprofile[i]);
  }
  for (i = 0; i < proot->height; i++)
  {
    print_next = 0;
    print_level(&print_next, proot, -xmin, i);
    printf("\n");
  }
  if (proot->height >= MAX_HEIGHT)
  {
    printf("(This tree is taller than %d, and may be drawn incorrectly.)\n",
           MAX_HEIGHT);
  }
  free_ascii_tree(proot);
}

//This function prints the given level of the given tree, assuming
//that the node has the given x cordinate.
void print_level(int *print_next, asciinode *node, int x, int level)
{
  int i, isleft;
  if (node == NULL) return;
  isleft = (node->parent_dir == -1);
  if (level == 0)
  {
    for (i=0; i<(x - (*print_next) -((node->lablen-isleft)/2)); i++)
    {
      printf(" ");
    }
    (*print_next) += i;
    printf("%s", node->label);
    (*print_next) += node->lablen;
  }
  else if (node->edge_length >= level)
  {
    if (node->left != NULL)
    {
      for (i=0; i<(x - (*print_next) - (level)); i++)
      {
        printf(" ");
      }
      (*print_next) += i;
      printf("/");
      (*print_next)++;
    }
    if (node->right != NULL)
    {
      for (i=0; i<(x - *print_next + (level)); i++)
      {
        printf(" ");
      }
      (*print_next) += i;
      printf("\\");
      (*print_next)++;
    }
  }
  else
  {
    print_level(print_next, node->left,
                    x-node->edge_length-1,
                    level-node->edge_length-1);
    print_level(print_next, node->right,
                    x+node->edge_length+1,
                    level-node->edge_length-1);
  }
}

int MIN (int X, int Y)
{
  return ((X) < (Y)) ? (X) : (Y);
}

int MAX (int X, int Y)
{
  return ((X) > (Y)) ? (X) : (Y);
}



//This function fills in the edge_length and
//height fields of the specified tree
void compute_edge_lengths(int gap, int *lprofile, int *rprofile, asciinode *node)
{
  int h, hmin, i, delta;
  if (node == NULL) return;
  compute_edge_lengths(gap, lprofile, rprofile, node->left);
  compute_edge_lengths(gap, lprofile, rprofile, node->right);

      /* first fill in the edge_length of node */
      if (node->right == NULL && node->left == NULL)
      {
        node->edge_length = 0;
      }
      else
      {
        if (node->left != NULL)
        {
          for (i=0; i<node->left->height && i < MAX_HEIGHT; i++)
          {
            rprofile[i] = -INFIN;
          }
          compute_rprofile(rprofile, node->left, 0, 0);
          hmin = node->left->height;
        }
        else
        {
          hmin = 0;
        }
        if (node->right != NULL)
        {
          for (i=0; i<node->right->height && i < MAX_HEIGHT; i++)
          {
            lprofile[i] = INFIN;
          }
          compute_lprofile(lprofile, node->right, 0, 0);
          hmin = MIN(node->right->height, hmin);
        }
        else
        {
          hmin = 0;
        }
        delta = 4;
        for (i=0; i<hmin; i++)
        {
          delta = MAX(delta, gap + 1 + rprofile[i] - lprofile[i]);
        }

        //If the node has two children of height 1, then we allow the
        //two leaves to be within 1, instead of 2
        if (((node->left != NULL && node->left->height == 1) ||
              (node->right != NULL && node->right->height == 1))&&delta>4)
        {
          delta--;
        }

        node->edge_length = ((delta+1)/2) - 1;
      }

      //now fill in the height of node
      h = 1;
      if (node->left != NULL) 
      {
        h = MAX(node->left->height + node->edge_length + 1, h);
      }
      if (node->right != NULL) 
      {
        h = MAX(node->right->height + node->edge_length + 1, h);
      }
      node->height = h;
    }

    asciinode * build_ascii_tree_recursive(IMG_NODE_t * t) 
    {
      asciinode * node;

      if (t == NULL) return NULL;

      node = malloc(sizeof(asciinode));
      node->left = build_ascii_tree_recursive(t->children);
      node->right = build_ascii_tree_recursive(t->sib);

      if (node->left != NULL)
      {
        node->left->parent_dir = -1;
      }

      if (node->right != NULL)
      {
        node->right->parent_dir = 1;
      }

      sprintf(node->label, "%d", t->to_image->frag.lb);
      node->lablen = strlen(node->label);

      return node;
    }


    //Copy the tree into the ascii node structre
    asciinode * build_ascii_tree(IMG_NODE_t * t) 
    {
      asciinode *node;
      if (t == NULL) return NULL;
      node = build_ascii_tree_recursive(t);
      node->parent_dir = 0;
      return node;
    }

    //Free all the nodes of the given tree
    void free_ascii_tree(asciinode *node) 
    {
      if (node == NULL) return;
      free_ascii_tree(node->left);
      free_ascii_tree(node->right);
      free(node);
    }

    //The following function fills in the lprofile array for the given tree.
    //It assumes that the center of the label of the root of this tree
    //is located at a position (x,y).  It assumes that the edge_length
    //fields have been computed for this tree.
    void compute_lprofile(int *lprofile, asciinode *node, int x, int y)
    {
      int i, isleft;
      if (node == NULL) return;
      isleft = (node->parent_dir == -1);
      lprofile[y] = MIN(lprofile[y], x-((node->lablen-isleft)/2));
      if (node->left != NULL) 
      {
        for (i=1; i <= node->edge_length && y+i < MAX_HEIGHT; i++) 
        {
          lprofile[y+i] = MIN(lprofile[y+i], x-i);
        }
      }
      compute_lprofile(lprofile, node->left, x-node->edge_length-1, y+node->edge_length+1);
      compute_lprofile(lprofile, node->right, x+node->edge_length+1, y+node->edge_length+1);
    }

    void compute_rprofile(int *rprofile, asciinode *node, int x, int y)
    {
      int i, notleft;
      if (node == NULL) return;
      notleft = (node->parent_dir != -1);
      rprofile[y] = MAX(rprofile[y], x+((node->lablen-notleft)/2));
      if (node->right != NULL) 
      {
        for (i=1; i <= node->edge_length && y+i < MAX_HEIGHT; i++) 
        {
          rprofile[y+i] = MAX(rprofile[y+i], x+i);
        }
      }
      compute_rprofile(rprofile, node->left, x-node->edge_length-1, y+node->edge_length+1);
      compute_rprofile(rprofile, node->right, x+node->edge_length+1, y+node->edge_length+1);
    }



