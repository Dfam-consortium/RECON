/*
 * ele_print.c  --  toString-style debug print helpers for ele data structures
 *
 * Implements all print_* functions whose prototypes are declared in ele.h.
 * These are used only for debugging and are never called in production paths.
 *
 * Author: Zhirong Bao
 * Modifications: Robert Hubley, Institute for Systems Biology
 */

#include "ele.h"
#include "ele_db.h"


/* ============================================================
 * Element, edge, and family print helpers
 * ============================================================ */

void print_ele_data(ELE_DATA_t *rt) {
  if (!rt) return;
  print_ele_info( rt->ele_info );
  if (rt->next)
    print_ele_data(rt->next);
}

void print_ele_info(ELE_INFO_t *rt) {
  if (!rt) return;
  printf("ELE_INFO_t: index=%d, stat=%c, file_update=%d,", rt->index, rt->stat,
         rt->file_updated );
  if ( rt->to_family )
    printf(" to_family->index=%d,", rt->to_family->index );
  else
    printf(" to_family=Null,");
  if ( rt->ele )
    printf(" ele->index=%d", rt->ele->index );
  else
    printf(" ele=Null");
  printf("\n");
}

void print_ele(ELEMENT_t *rt) {
  if (!rt) return;
  printf("ELEMENT_t: index=%d, direction=%d, update=%d, l_hold=%d, img_no=%d, flimg_no=%d, edge_no=%d,",
         rt->index, rt->direction, rt->update, rt->l_hold, rt->img_no, rt->flimg_no, rt->edge_no );
  if ( rt->to_img_tree )
    printf(" to_img_tree->to_image->index=%d,", rt->to_img_tree->to_image->index );
  else
    printf(" to_img_tree=Null,");
  if ( rt->to_img_data )
    printf(" to_img_data->to_image->index=%d,", rt->to_img_data->to_image->index );
  else
    printf(" to_img_data=Null,");
  if ( rt->edges )
    printf(" edges->to_edge->index=%d,", rt->edges->to_edge->index );
  else
    printf(" edges=Null,");
  if ( rt->PCP )
    printf(" PCP->cp=%d,", rt->PCP->cp );
  else
    printf(" PCP=Null,");
  if ( rt->TBD )
    printf(" TBD->bd=%d,", rt->TBD->bd );
  else
    printf(" TBD=Null,");
  if ( rt->redef )
    printf(" redef->ele_info->index=%d", rt->redef->ele_info->index );
  else
    printf(" redef=Null");
  printf("\n");
}

void print_edge(EDGE_t *rt) {
  if (!rt) return;
  printf("EDGE_t: index=%d, type=%c, direction=%d, score=%d,", rt->index, rt->type,
         rt->direction, rt->score );
  if ( rt->ele1_info )
    printf(" ele1_info->index=%d,", rt->ele1_info->index );
  else
    printf(" ele1_info=Null,");
  if ( rt->ele2_info )
    printf(" ele2_info->index=%d,", rt->ele2_info->index );
  else
    printf(" ele2_info=Null,");
  printf("\n");
}

/*
 * print_edge_tree_GML  --  DFS of an edge tree, printing edges in GML format.
 *
 * If rel_to_ele_id < 1, print all edges in the tree.
 * Otherwise only print edges where the lower-index endpoint matches
 * rel_to_ele_id (produces a non-redundant set of edges).
 */
void print_edge_tree_GML(EDGE_TREE_t *rt, int rel_to_ele_id) {
  if (!rt) return;
  if ( rel_to_ele_id < 1 ||
       ((rt->to_edge->ele1_info->index < rt->to_edge->ele2_info->index &&
        rt->to_edge->ele1_info->index == rel_to_ele_id ) ||
       (rt->to_edge->ele2_info->index < rt->to_edge->ele1_info->index &&
        rt->to_edge->ele2_info->index == rel_to_ele_id )) )
    {
    printf("edge [\n");
    printf("       source %d\n", rt->to_edge->ele1_info->index );
    printf("       target %d\n", rt->to_edge->ele2_info->index );
    printf("       label \"index=%d, type=%c, direction=%d, score=%d\"\n", rt->to_edge->index, rt->to_edge->type,
         rt->to_edge->direction, rt->to_edge->score );
    printf("     ]\n");
  }
  if (rt->l)
    print_edge_tree_GML(rt->l, rel_to_ele_id);
  if (rt->r)
    print_edge_tree_GML(rt->r, rel_to_ele_id);
}


void print_all_eles_GML() {
  int i = 0;
  printf("graph [\n");
  for (i=0; i<ele_ct && i<ele_array_size; i++) {
    printf("node [\n");
    printf("       id %d\n", (*(all_ele+i))->index);
    printf("       label \"e%d\"\n", (*(all_ele+i))->index);
    printf("     ]\n");
  }
  for (i=0; i<ele_ct && i<ele_array_size; i++)
    print_edge_tree_GML((*(all_ele+i))->ele->edges, (*(all_ele+i))->ele->index);
  printf("]\n");
}


void print_local_network(ELE_DATA_t *rt) {
  if (!rt) return;
  ELE_DATA_t *curr;
  curr=rt;
  while ( curr ) {
    print_edge_tree_GML(curr->ele_info->ele->edges, curr->ele_info->ele->index);
    curr = curr->next;
  }
}

void print_edge_tree(EDGE_TREE_t *rt, int level) {
  int i;
  if (!rt) return;
  for ( i = 0; i < level; i++ )
    printf("  ");
  print_edge(rt->to_edge);
  if (rt->l)
    print_edge_tree(rt->l, level+1);
  if (rt->r)
    print_edge_tree(rt->r, level+1);
}


/* ============================================================
 * CP, BD, and family print helpers
 * ============================================================ */

void print_cp(CP_t *cp) {
  if (!cp) { printf("CP_t: NULL\n"); return; }
  printf("CP_t: cp=%d, contributor->index=%d (stat=%c)\n",
         cp->cp,
         cp->contributor ? cp->contributor->index : -1,
         cp->contributor ? cp->contributor->stat  : '?');
}

void print_cp_list(CP_t *cp) {
  int count = 0;
  while (cp) {
    printf("  [%d] ", count++);
    print_cp(cp);
    cp = cp->next;
  }
}

void print_bd(BD_t *bd) {
  if (!bd) { printf("BD_t: NULL\n"); return; }
  printf("BD_t: bd=%d, support=%d\n", bd->bd, bd->support);
}

void print_bd_list(BD_t *bd) {
  int count = 0;
  while (bd) {
    printf("  [%d] ", count++);
    print_bd(bd);
    bd = bd->next;
  }
}

void print_family(FAMILY_t *fam) {
  ELE_DATA_t *m;
  int mem_count = 0;
  if (!fam) { printf("FAMILY_t: NULL\n"); return; }
  m = fam->members;
  while (m) { mem_count++; m = m->next; }
  printf("FAMILY_t: index=%d, name='%s', member_count=%d\n",
         fam->index, fam->name, mem_count);
}

void print_fam_data(FAM_DATA_t *fd) {
  int count = 0;
  while (fd) {
    printf("  [%d] ", count++);
    print_family(fd->to_family);
    fd = fd->next;
  }
}
