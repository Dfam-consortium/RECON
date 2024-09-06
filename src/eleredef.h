#ifndef __ELEREDEF_H__
#define __ELEREDEF_H__


typedef struct img_node {
  short recorded;
  IMAGE_t *to_image;
  struct img_node *sib;
  struct img_node *children;
} IMG_NODE_t;

// Function signatures
int print_consis_tree_GV( IMG_NODE_t *rt, int parent_index, int index );
int print_consis_tree(IMG_NODE_t *rt);


// Debug function to view consistency tree in GraphViz format
int print_consis_tree_GV( IMG_NODE_t *rt, int parent_index, int index ) {
  if (!rt) return index;

  if ( index == 0 ) {
    // Top level:
    parent_index = 0;
    printf("CONSISTENCY TREE GV-format:\n");

    if ( rt->sib != NULL )
      printf(" -- warning: top-level token node has siblings!\n");

    printf("digraph consis_tree {\n");
    printf("  graph [ label=\"Consistency Tree for E%d to E%d\" ]\n", rt->children->to_image->ele_info->index, partner(rt->children->to_image)->ele_info->index);
    printf("    idx0 [label=\"root\"]\n");
    printf("    idx0\n");
    index = print_consis_tree_GV( rt->children, parent_index, index + 1 );
    printf("}\n");
    return index;
  }

  IMAGE_t *ipt = partner(rt->to_image);
  int node_index = index;
  if ( rt->to_image->to_msp->direction == 1 )
    printf("    idx%d [label=\"%d-%d/%d-%d[+]\"", index, rt->to_image->frag.lb, rt->to_image->frag.rb, ipt->frag.lb, ipt->frag.rb );
  else
    printf("    idx%d [label=\"%d-%d/%d-%d[-]\"", index, rt->to_image->frag.lb, rt->to_image->frag.rb, ipt->frag.lb, ipt->frag.rb );

  if ( rt->to_image->to_msp->stat == 'p' )
    printf(", style=filled, fillcolor=lightblue]\n");
  else
    printf("]\n");

  //if ( rt->sib != NULL && rt->children != NULL )
  //  printf("    idx%d [label=\"p%d-%d/%d-%d[c,s]\"]\n", index, rt->to_image->frag.lb, rt->to_image->frag.rb, ipt->frag.lb, ipt->frag.rb );
  //else if ( rt->sib != NULL && rt->children == NULL )
  //  printf("    idx%d [label=\"p%d-%d/%d-%d[s]\"]\n", index, rt->to_image->frag.lb, rt->to_image->frag.rb, ipt->frag.lb, ipt->frag.rb );
  //else if ( rt->sib == NULL && rt->children != NULL )
  //  printf("    idx%d [label=\"p%d-%d/%d-%d[c]\"]\n", index, rt->to_image->frag.lb, rt->to_image->frag.rb, ipt->frag.lb, ipt->frag.rb );
  //else
  //  printf("    idx%d [label=\"p%d-%d/%d-%d[]\"]\n", index, rt->to_image->frag.lb, rt->to_image->frag.rb, ipt->frag.lb, ipt->frag.rb );

  printf("    idx%d -> idx%d\n", parent_index, index );
  index++;

  if ( rt->sib ) {
    index = print_consis_tree_GV( rt->sib, parent_index, index );
  }

  if ( rt->children != NULL ){
    index = print_consis_tree_GV( rt->children, node_index, index );
  }

  return index;

}

// RMH - debug function to view consistency tree
int print_consis_tree(IMG_NODE_t *rt) {
  int level = 0;
  IMG_NODE_t *child_rt;
  IMG_NODE_t *sib_rt;

  // Top level:
  printf("CONSISTENCY TREE:\n");
  if ( rt->sib != NULL )
    printf("  c_tree: level=%d [(null):0:0-0] ( has sibs - unexpected! )\n",level);
  else
    printf("  c_tree: level=%d [(null):0:0-0]\n",level);

  level++;
  child_rt = rt->children;
  while ( child_rt != NULL ) {
    IMAGE_t *ipt = partner(child_rt->to_image);
    printf("  c_tree: level=%d e%d:%s:%d:%d-%d to e%d:%s:%d:%d-%d : siblings[",level, child_rt->to_image->ele_info->index, child_rt->to_image->frag.seq_name, child_rt->to_image->to_msp->direction,  child_rt->to_image->frag.lb, child_rt->to_image->frag.rb, ipt->ele_info->index, ipt->frag.seq_name, ipt->to_msp->direction, ipt->frag.lb, ipt->frag.rb);
    sib_rt = child_rt->sib;
    while ( sib_rt != NULL )
    {
      IMAGE_t *ipt = partner(sib_rt->to_image);
      if ( sib_rt->children == NULL )
        printf("e%d:%s:%d:%d-%d to e%d:%s:%d:%d-%d,", sib_rt->to_image->ele_info->index, sib_rt->to_image->frag.seq_name, sib_rt->to_image->to_msp->direction, sib_rt->to_image->frag.lb, sib_rt->to_image->frag.rb, ipt->ele_info->index, ipt->frag.seq_name, ipt->to_msp->direction, ipt->frag.lb, ipt->frag.rb);
      else
        printf("%s:%d:%d-%d(*e%d:%d),", sib_rt->to_image->frag.seq_name, sib_rt->to_image->to_msp->direction, sib_rt->to_image->frag.lb, sib_rt->to_image->frag.rb,
                                        sib_rt->children->to_image->ele_info->index, sib_rt->children->to_image->frag.lb);
      sib_rt = sib_rt->sib;
    }
    printf("]\n");
    level++;
    child_rt = child_rt->children;
  }
  return 0;
}


#endif
