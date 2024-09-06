/*
 * eledef.c
 *
 * Element definition using single linkage clustering
 *
 * Author: Zhirong Bao
 * Minor modifications by: Robert Hubley, Institute for Systems Biology
 *
 * Description:
 * This program reads sequence and MSP data to perform clustering based on single or double linkage
 * methods. It processes command-line inputs, handles sequence names and MSPs, and outputs results
 * to various files. Unit tests are provided using the minunit framework.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <stdint.h>
#include <getopt.h>
#include <minunit.h>

// Constants
#define NAME_LEN 50
#define MIN_OVERLAP 10
#define IMG_CAP 500000

// Struct definitions
typedef struct frag {  // msps.h ( modified to uint64_t )
    char *seq_name;
    uint64_t lb, rb;
} FRAG_t;

typedef struct frag_list { // msps.h
    FRAG_t *to_frag;
    struct frag_list *next;
} FRAG_DATA_t;

typedef struct image { // msps.h
    int index;
    FRAG_t frag;
    struct msp *to_msp;
    struct ele_info *ele_info;
} IMAGE_t;

typedef struct image_list { // msps.h
    IMAGE_t *to_image;
    struct image_list *next;
} IMG_DATA_t;

typedef struct msp { // msps.h
    char stat; // 'p' = highest scoring MSP for a given element pair it must be full length
               //       for at least one of the elements.
               // 's'
    uint32_t score;
    float iden;
    int direction;
    IMAGE_t query, sbjct;
} MSP_t;

typedef struct msp_list { // msps.h
    MSP_t *to_msp;
    struct msp_list *next;
} MSP_DATA_t;

typedef struct msp_prototype {
    int pe, se;
} MPROT_t;

typedef struct img_prototype {
    int index, ele_index;
    MSP_t *to_msp;
} IPROT_t;

typedef struct ele_prototype {
    short flag; // 0 = not printed, 1 = printed
    int index, img_no;
    FRAG_t frag;
    struct ele_prototype *next;
} EPROT_t;

// Function prototypes
int get_seq_names(FILE *seq_list, char ***seq_names);
int get_seq_index(int, int, char *, char **);
int print_msp(MSP_t *);
int fprint_msp(FILE *, MSP_t *);
int scan_msp(MSP_t *, char *, char **, int);
IMAGE_t *partner(IMAGE_t *);
int doub_cov(FRAG_t *, FRAG_t *, float);
int sing_cov(FRAG_t *, FRAG_t *, float);
void ele_def(int, FILE *, FILE *, float, EPROT_t **, int *, MPROT_t **, FILE *, FILE *, char **, int);
void img_charge(IPROT_t **, int, FILE *, char **, int);
int index_cmp(const void *, const void *);
void include_image(int, IMAGE_t *, FRAG_t *, int *, float, FILE *, MPROT_t **, IMG_DATA_t **, IMG_DATA_t **, int *);
void save_element_to_ep_list(EPROT_t **, EPROT_t **, int *, int, FRAG_t, FILE *);
void init_element_from_remain_list(int *, int *, FRAG_t *, IMAGE_t *, IMG_DATA_t **, IMG_DATA_t **, IMG_DATA_t **, FILE *, MPROT_t **, float);
void save_img_to_remaining_list(IMAGE_t, IMG_DATA_t **, IMG_DATA_t **);
void build_element_files(int, IPROT_t **, IPROT_t **, FILE *, MPROT_t **, EPROT_t **, FILE **, FILE *, FILE *, char **, int);


/**
 * @brief Reads sequence names from a file and returns them along with the count.
 *
 * This function reads a list of sequence names from a given file and stores them
 * in a dynamically allocated array. It also returns the count of sequences read.
 *
 * This file *MUST* be in lexically sorted order as per the RECON documentation.
 * The MSP file must also be in lexically sorted order (by query) sequence.
 * TODO: Remove this requirement by sorting at the start of the program.
 *
 * @param seq_list File pointer to the sequence list.
 * @param seq_names Pointer to an array where sequence names will be stored.
 * @return The number of sequences read. Returns -1 if there's an error.
 */
int get_seq_names(FILE *seq_list, char ***seq_names) {
    char line[256];
    int seq_ct = 0;
    char *name_start;
    int seq_no;

    // Get the number of sequences in the list, and generates an array of proper size to hold the list
    if (fgets(line, sizeof(line), seq_list)) {
        seq_no = atoi(line);
        if (seq_no <= 0) {
            fprintf(stderr, "First line in the list of sequences must be a positive integer. Exit\n");
            return -1;
        }
    } else {
        fprintf(stderr, "Failed to read from file. Exit\n");
        return -1;
    }

    // Allocate memory for sequence names
    *seq_names = (char **)malloc(seq_no * sizeof(char *));
    if (*seq_names == NULL) {
        fprintf(stderr, "Memory allocation failed. Exit\n");
        return -1;
    }

    // Read sequence names
    while (fgets(line, sizeof(line), seq_list) && seq_ct < seq_no) {
        (*seq_names)[seq_ct] = (char *)malloc(NAME_LEN * sizeof(char));
        if ((*seq_names)[seq_ct] == NULL) {
            fprintf(stderr, "Memory allocation failed. Exit\n");
            // Free previously allocated memory
            for (int i = 0; i < seq_ct; i++) {
                free((*seq_names)[i]);
            }
            free(*seq_names);
            return -1;
        }

        // Skip '>' and any following spaces if present
        name_start = line;
        if (line[0] == '>') {
            name_start++;
            while (isspace((unsigned char)*name_start)) name_start++;
        }

        // Copy name up to NAME_LEN-1 characters
        strncpy((*seq_names)[seq_ct], name_start, NAME_LEN - 1);
        (*seq_names)[seq_ct][NAME_LEN - 1] = '\0';  // Ensure null-terminated string

        // Trim trailing whitespace
        name_start = (*seq_names)[seq_ct];
        while (*name_start && !isspace((unsigned char)*name_start)) name_start++;
        *name_start = '\0';

        seq_ct++;
    }

    return seq_ct;
}

/**
 * @brief Binary search of a sorted list of sequence names.
 *
 * @param left Left index.
 * @param right Right index.
 * @param seq_name Sequence name to search for.
 * @param seq_names A **lexically sorted** array of sequence names.
 * @return Position in the array if found, -1 if not found.
 */
int get_seq_index(int left, int right, char *seq_name, char **seq_names) {
    int pos, dir;

    if (left <= right) {
        pos = (left + right) / 2;
        dir = strncmp(seq_name, seq_names[pos], NAME_LEN);
        if (dir < 0) return get_seq_index(left, pos - 1, seq_name, seq_names);
        else if (dir > 0) return get_seq_index(pos + 1, right, seq_name, seq_names);
        else return pos;
    } else {
        return -1;
    }
}

/**
 * @brief Prints an MSP to stdout.
 * 
 * @param m MSP structure.
 * @return 0 if successful, 1 otherwise.
 */
int print_msp(MSP_t *m) {
    if (m->direction == 1) {
        return printf("%06d %3.1f %010ld %010ld %s %010ld %010ld %s \n", m->score, m->iden,
                      m->query.frag.lb, m->query.frag.rb, m->query.frag.seq_name,
                      m->sbjct.frag.lb, m->sbjct.frag.rb, m->sbjct.frag.seq_name) > 0 ? 0 : 1;
    } else if (m->direction == -1) {
        return printf("%06d %3.1f %010ld %010ld %s %010ld %010ld %s \n", m->score, m->iden,
                      m->query.frag.rb, m->query.frag.lb, m->query.frag.seq_name,
                      m->sbjct.frag.lb, m->sbjct.frag.rb, m->sbjct.frag.seq_name) > 0 ? 0 : 1;
    } else {
        return 1;
    }
}

/**
 * @brief Writes an MSP to a file.
 * 
 * @param ofp Output file pointer.
 * @param m MSP structure.
 * @return 0 if successful, 1 otherwise.
 */
int fprint_msp(FILE *ofp, MSP_t *m) {
    if (m->direction == 1) {
        return fprintf(ofp, "%06d %3.1f %010ld %010ld %s %010ld %010ld %s \n", m->score, m->iden,
                       m->query.frag.lb, m->query.frag.rb, m->query.frag.seq_name,
                       m->sbjct.frag.lb, m->sbjct.frag.rb, m->sbjct.frag.seq_name) > 0 ? 0 : 1;
    } else if (m->direction == -1) {
        return fprintf(ofp, "%06d %3.1f %010ld %010ld %s %010ld %010ld %s \n", m->score, m->iden,
                       m->query.frag.rb, m->query.frag.lb, m->query.frag.seq_name,
                       m->sbjct.frag.lb, m->sbjct.frag.rb, m->sbjct.frag.seq_name) > 0 ? 0 : 1;
    } else {
        return 1;
    }
}

/**
 * @brief Scans an MSP line and initializes an MSP structure.
 *
 *  -- Also in msps.h
 *
 * @param m MSP structure to initialize.
 * @param line Line containing MSP data.
 * @param seq_names Array of sequence names.
 * @param seq_no Number of sequences.
 * @return 0 if successful, 1 otherwise.
 */
int scan_msp(MSP_t *m, char *line, char **seq_names, int seq_no) {
    int32_t bd_tmp;
    char qname[NAME_LEN], sname[NAME_LEN];
    int pos;

    if (sscanf(line, "%d %f %ld %ld %s %ld %ld %s \n", &(m->score), &(m->iden),
               &(m->query.frag.lb), &(m->query.frag.rb), qname,
               &(m->sbjct.frag.lb), &(m->sbjct.frag.rb), sname) != 8) {
        return 1;
    } else {
        pos = get_seq_index(0, seq_no - 1, qname, seq_names);
        if (pos < 0) return pos;
        m->query.frag.seq_name = seq_names[pos];
        pos = get_seq_index(0, seq_no - 1, sname, seq_names);
        if (pos < 0) return pos;
        m->sbjct.frag.seq_name = seq_names[pos];

        m->direction = 1;
        m->stat = 's';
        m->query.to_msp = m;
        m->sbjct.to_msp = m;

        if (m->query.frag.lb > m->query.frag.rb) {
            m->direction *= -1;
            bd_tmp = m->query.frag.lb;
            m->query.frag.lb = m->query.frag.rb;
            m->query.frag.rb = bd_tmp;
        }
        if (m->sbjct.frag.lb > m->sbjct.frag.rb) {
            m->direction *= -1;
            bd_tmp = m->sbjct.frag.lb;
            m->sbjct.frag.lb = m->sbjct.frag.rb;
            m->sbjct.frag.rb = bd_tmp;
        }
        return 0;
    }
}

/**
 * @brief Returns the partner image in an MSP.
 *
 * @param i IMAGE structure.
 * @return Partner IMAGE structure.
 */
IMAGE_t *partner(IMAGE_t *i) {
    return (i == &(i->to_msp->query)) ? &(i->to_msp->sbjct) : &(i->to_msp->query);
}

/**
 * @brief Check single coverage between two fragments.
 *
 * Criteria:
 *     o Same sequence name.
 *     o Coverage of either fragment is greater than the cutoff
 *       fraction: cutoff * length of the fragment.
 *
 *  -- Also in msps.h
 *
 * @param f1 First fragment.
 * @param f2 Second fragment.
 * @param cutoff Cutoff for coverage.
 * @return 1 if coverage meets the criteria, 0 otherwise.
 */
int sing_cov(FRAG_t *f1, FRAG_t *f2, float cutoff) {
    //int32_t l1, l2, l, lb, rb;
    uint64_t l1, l2, l, lb, rb;
    if (f1->seq_name == f2->seq_name) {
        l1 = f1->rb - f1->lb;
        l2 = f2->rb - f2->lb;
        lb = f1->lb > f2->lb ? f1->lb : f2->lb;
        rb = f1->rb < f2->rb ? f1->rb : f2->rb;
        l = rb - lb;
        if ((float) l / l1 >= cutoff || (float) l / l2 >= cutoff) {
            return 1;
        }
    }
    return 0;
}

/**
 * @brief Check double coverage between two fragments.
 *
 *  -- Also in msps.h
 *
 * @param f1 First fragment.
 * @param f2 Second fragment.
 * @param cutoff Cutoff for coverage.
 * @return 1 if coverage meets the criteria, 0 otherwise.
 */
int doub_cov(FRAG_t *f1, FRAG_t *f2, float cutoff) {
    //int32_t l1, l2, l, lb, rb;
    uint64_t l1, l2, l, lb, rb;
    if (f1->seq_name == f2->seq_name) {
        l1 = f1->rb - f1->lb;
        l2 = f2->rb - f2->lb;
        lb = f1->lb > f2->lb ? f1->lb : f2->lb;
        rb = f1->rb < f2->rb ? f1->rb : f2->rb;
        l = rb - lb;
        if ((float) l / l1 >= cutoff && (float) l / l2 >= cutoff) {
            return 1;
        }
    }
    return 0;
}

/**
 * @brief Load MSP referenced by image structure entry
 *
 * @param shadow Shadow array of IPROT_t pointers.
 * @param ct Count of images.
 * @param input Input file pointer.
 */
void img_charge(IPROT_t **shadow, int ct, FILE *input, char **seq_names, int seq_no) {
    int i = 0, pos = 0;
    char line[151];
    int scan_flag;
    MSP_t msp;

    qsort(shadow, ct, sizeof(IPROT_t *), index_cmp);
    while (fgets(line, 150, input)) {
        if (pos == (*(shadow + i))->index / 2) {
            scan_flag = scan_msp(&msp, line, seq_names, seq_no);
            if (scan_flag) {
                (*(shadow + i))->to_msp->score = 0;
                fprintf(stderr, "Wrong format in the MSP file line %d for image %d\n", pos, (*(shadow + i))->index);
            } else {
                *((*(shadow + i))->to_msp) = msp;
            }
            i++;
            if (i == ct) break;
            if (pos == (*(shadow + i))->index / 2) {
                // In case the next image is from the same MSP line
                if (scan_flag) {
                    (*(shadow + i))->to_msp->score = 0;
                    fprintf(stderr, "Wrong format in the MSP file line %d for image %d\n", pos, (*(shadow + i))->index);
                } else {
                    *((*(shadow + i))->to_msp) = msp;
                }
                i++;
                if (i == ct) break;
            }
        }
        pos++;
    }
    rewind(input);
}

/**
 * @brief Compare indices of IPROT_t for sorting.
 *
 * @param i1 First index.
 * @param i2 Second index.
 * @return Difference between indices.
 */
int index_cmp(const void *i1, const void *i2) {
    return (*((IPROT_t **) i1))->index - (*((IPROT_t **) i2))->index;
}

/**
 * @brief Defines elements using linkage clustering.
 *
 * @param method Clustering method.
 * @param frags Input file pointer for fragments.
 * @param img_prot file for saving images.
 * @param cutoff Cutoff value for clustering.
 * @param all_epp Pointer to list of all elements.
 * @param ecp Pointer to element count.
 * @param all_mprot Array of all MSP prototypes.
 * @param all_ele Output file pointer for all elements.
 * @param ele_no Output file pointer for element count.
 */
void ele_def(int method, FILE *frags, FILE *img_prot, float cutoff, EPROT_t **all_epp, int *ecp, MPROT_t **all_mprot, FILE *all_ele, FILE *ele_no, char **seq_names, int seq_no) {
    int i, img_ct, cov;
    char line[100];
    FRAG_t ele_frag;
    char fragname[NAME_LEN];
    int pos;
    IMAGE_t img, *imgp;
    IMG_DATA_t *cur, *prev = NULL, *remain = NULL, *tail = NULL, *tmp;
    int ritetime;
    EPROT_t *ep_tail, *ep_tmp;

    // ritetime marks the right time to start a new element
    ritetime = 1;

    // Read in select columns from the images_sorted file:
    //   - Image Index/img.index: This is a unique identifier for the image and is translatable
    //                            back to the source MSP line
    //   - Sequence Identifier/fragname: The unique identifier for the IMAGE sequence.
    //   - Left Boundary/frag.lb: The start position of the IMAGE in the
    //                            sequence (start<end, 1 based, fully closed)
    //   - Right Boundary/frag.rb: The end position of the IMAGE in the
    //                             sequence (start<end, 1 based, fully closed)
    //
    //  Algorithm:
    //      Loop over images
    //          - Initialize a new element given current image if previous element complete (ritetime)
    //            then continue to next image.
    //          - Check for image overlap with current element
    //          - If overlap by MIN_OVERLAP:
    //                - If meets single linkage clustering criteria:
    //                   - Add image to element
    //                    - if image is longer than current element:
    //                        - Update element right boundary
    //                        - Check remaining linked list for disqualified images to include with
    //                          updated boundary.
    //                  Else (single linage fail):
    //                    - Save image to remaining linked list
    //            Else (<minimum overlap)
    //                - Finalize element
    //                - Save current image to remaining linked list
    //                - Set element complete flag (ritetime)
    //                - Loop over remaining linked list until it's used up and then continue to
    //                  next image in the main loop.  ( see note about this complexity not being
    //                  necessary below ).
    //
    while (fgets(line, 100, frags)) {
        sscanf(line, "%d %*d %s %ld %ld\n", &img.index, fragname, &img.frag.lb, &img.frag.rb);
        pos = get_seq_index(0, seq_no - 1, fragname, seq_names);
        img.frag.seq_name = seq_names[pos];

        if (ritetime) {
            ritetime = 0;
            // Increase the element count (through the passed pointer)
            (*ecp)++;
            img_ct = 1;
            // Start the definition of the element
            ele_frag = img.frag;

            // Fill in the map from MSP index to element....note 0=no element mapping, and 1..n are element
            // numbers.
            if (img.index % 2) {
                // Subject range of the MSP
                (*(all_mprot + img.index / 2))->se = (*ecp);
            } else {
                // Query range of the MSP
                (*(all_mprot + img.index / 2))->pe = (*ecp);
            }
            // Save Element#, Image#
            fprintf(img_prot, "%d %d\n", (*ecp), img.index);
            fflush(img_prot);
            continue;
        }

        // Checking if sequences overlap by more than MIN_OVERLAP
        //
        //       lb|--ele_frag--|rb
        //                  lb|--img.frag--|rb
        //
        //   ( NOTE: This is a pointer comparison rather than
        //   a string comparison -- be safe out there )
        //  RMH: Original code had this minimum overlap hardcoded
        //       to 10bp.
        if (ele_frag.seq_name == img.frag.seq_name && ele_frag.rb - img.frag.lb > MIN_OVERLAP) {
            // RMH: sing_cov from msps.h
            //      Return 1 if the same sequence *and* the overlap
            //             between the two sequences is >= cutoff * length
            //             of either sequence.
            cov = (method == 1) ? sing_cov(&ele_frag, &img.frag, cutoff) : doub_cov(&ele_frag, &img.frag, cutoff);
            if (cov) {
                // A good image; include it in the current element
                img_ct++;
                if (img.index % 2) {
                    // Odd indices originate from the subject range of an MSP
                    (*(all_mprot + img.index / 2))->se = (*ecp);
                } else {
                    // Even indices originate from the query range of an MSP
                    (*(all_mprot + img.index / 2))->pe = (*ecp);
                }
                fprintf(img_prot, "%d %d\n", (*ecp), img.index);
                fflush(img_prot);
                // If we find an overlapping image which is longer than the current
                // element range, update the element endpoint and reconsider
                // if previous unqualified images are now able to participate
                // in this element.
                if (ele_frag.rb < img.frag.rb) {
                    ele_frag.rb = img.frag.rb;
                    cur = remain;
                    prev = NULL;
                    while (cur) {
                        include_image(method, cur->to_image, &ele_frag, &img_ct, cutoff, img_prot, all_mprot, &remain, &tail, ecp);
                    }
                }
            } else {
                // Current image is not good for the current ele, keep it in the remaining
                // list and move on
                save_img_to_remaining_list(img, &remain, &tail);
            }
        } else {
            // Output the element
            save_element_to_ep_list(all_epp, &ep_tail, ecp, img_ct, ele_frag, all_ele);
            save_img_to_remaining_list(img, &remain, &tail);
            // Start a new element starting from the remain list
            ritetime = 1;

            // RMH:
            // I am not sure this is necessary.  I suspect this could would be more efficient
            // if this was redesigned as a buffered linked list of images.  The main loop would
            // simply always act on the start of the list and it would only get expanded as needed.
            while (remain && ritetime) {
                cur = remain;
                ritetime = 0;
                init_element_from_remain_list(ecp, &img_ct, &ele_frag, &img, &cur, &prev, &remain, img_prot, all_mprot, cutoff);
                // cur points to remain list at this point.  If
                // there are still more "remaining" try to include
                // them first.
                while (cur) {
                    if (ele_frag.seq_name == (*cur).to_image->frag.seq_name && ele_frag.rb - (*cur).to_image->frag.lb > MIN_OVERLAP) {
                        include_image(method, (*cur).to_image, &ele_frag, &img_ct, cutoff, img_prot, all_mprot, &remain, &tail, ecp);
                    } else {
                        save_element_to_ep_list(all_epp, &ep_tail, ecp, img_ct, ele_frag, all_ele);
                        ritetime = 1;
                        break;
                    }
                }
            }
        }
    }

    save_element_to_ep_list(all_epp, &ep_tail, ecp, img_ct, ele_frag, all_ele);
    if (remain) ritetime = 1;

    // Add the last element to the list of elements */
    while (remain) {
        cur = remain;
        prev = NULL;
        while (cur) {
            if (ritetime) {
                ritetime = 0;
                init_element_from_remain_list(ecp, &img_ct, &ele_frag, &img, &cur, &prev, &remain, img_prot, all_mprot, cutoff);
                continue;
            }
            if (ele_frag.seq_name == (*cur).to_image->frag.seq_name && ele_frag.rb - (*cur).to_image->frag.lb > MIN_OVERLAP) {
                include_image(method, (*cur).to_image, &ele_frag, &img_ct, cutoff, img_prot, all_mprot, &remain, &tail, ecp);
            } else {
                save_element_to_ep_list(all_epp, &ep_tail, ecp, img_ct, ele_frag, all_ele);
                ritetime = 1;
                break;
            }
        }
        if (!cur) {
            save_element_to_ep_list(all_epp, &ep_tail, ecp, img_ct, ele_frag, all_ele);
            ritetime = 1;
        }
    }

    fprintf(ele_no, "%d\n", (*ecp));
}

/**
 * @brief Includes an image in the element if it meets criteria.
 *
 * Re-evaluate orphan images (remains list) after an element boundary has been extended.
 *
 * @param method Clustering method.
 * @param image IMAGE structure to be included.
 * @param ele_frag Current element fragment.
 * @param img_ct Pointer to image count.
 * @param cutoff Cutoff value for clustering.
 * @param img_prot Image prototype file.
 * @param all_mprot Array of all MSP prototypes.
 * @param remain Pointer to the remaining list.
 * @param tail Pointer to the tail of the remaining list.
 * @param ecp Pointer to element count.
 */
void include_image(int method, IMAGE_t *image, FRAG_t *ele_frag, int *img_ct, float cutoff, FILE *img_prot, MPROT_t **all_mprot, IMG_DATA_t **remain, IMG_DATA_t **tail, int *ecp) {
    int cov = (method == 1) ? sing_cov(ele_frag, &image->frag, cutoff) : doub_cov(ele_frag, &image->frag, cutoff);
    if (cov) {
        (*img_ct)++;
        if (image->index % 2) {
            (*(all_mprot + image->index / 2))->se = (*ecp);
        } else {
            (*(all_mprot + image->index / 2))->pe = (*ecp);
        }
        fprintf(img_prot, "%d %d\n", (*ecp), image->index);
        fflush(img_prot);
        if (ele_frag->rb < image->frag.rb) {
            ele_frag->rb = image->frag.rb;
            free(image);
            *remain = NULL;
            *tail = NULL;
        } else {
            free(image);
            *remain = (*remain)->next;
            if (*remain) {
                *tail = *remain;
            }
        }
    } else {
        if (*tail) {
            (*tail)->next = *remain;
        } else {
            *remain = *remain;
        }
        *tail = *remain;
    }
}

/**
 * @brief Saves the current element to the EP list.
 *
 * @param all_epp Pointer to list of all elements.
 * @param ep_tail Pointer to the tail of the element list.
 * @param ecp Pointer to element count.
 * @param img_ct Image count.
 * @param ele_frag Element fragment.
 * @param all_ele All element file.
 */
void save_element_to_ep_list(EPROT_t **all_epp, EPROT_t **ep_tail, int *ecp, int img_ct, FRAG_t ele_frag, FILE *all_ele) {
    EPROT_t *ep_tmp = (EPROT_t *) malloc(sizeof(EPROT_t));
    ep_tmp->flag = 0;
    ep_tmp->index = (*ecp);
    ep_tmp->img_no = img_ct;
    ep_tmp->frag = ele_frag;
    ep_tmp->next = NULL;
    if (*all_epp) {
        (*ep_tail)->next = ep_tmp;
    } else {
        *all_epp = ep_tmp;
    }
    *ep_tail = ep_tmp;
    fprintf(all_ele, "%d %s %ld %ld\n", (*ecp), ele_frag.seq_name, ele_frag.lb, ele_frag.rb);
}

/**
 * @brief Initializes a new element from the remaining list.
 *
 * @param ecp Pointer to element count.
 * @param img_ct Pointer to image count.
 * @param ele_frag Pointer to element fragment.
 * @param img Pointer to image.
 * @param cur Pointer to the current list node.
 * @param prev Pointer to the previous list node.
 * @param remain Pointer to the remaining list.
 * @param img_prot Image prototype file.
 * @param all_mprot Array of all MSP prototypes.
 * @param cutoff Cutoff value for clustering.
 */
void init_element_from_remain_list(int *ecp, int *img_ct, FRAG_t *ele_frag, IMAGE_t *img, IMG_DATA_t **cur, IMG_DATA_t **prev, IMG_DATA_t **remain, FILE *img_prot, MPROT_t **all_mprot, float cutoff) {
    (*ecp)++;
    *img_ct = 1;
    *ele_frag = (*cur)->to_image->frag;

    if ((*cur)->to_image->index % 2) {
        (*(all_mprot + (*cur)->to_image->index / 2))->se = (*ecp);
    } else {
        (*(all_mprot + (*cur)->to_image->index / 2))->pe = (*ecp);
    }
    fprintf(img_prot, "%d %d\n", (*ecp), (*cur)->to_image->index);
    fflush(img_prot);
    *remain = (*cur)->next;
    free((*cur)->to_image);
    free(*cur);
    *prev = NULL;
    *cur = *remain;
}

/**
 * @brief Saves an image to the remaining list.
 *
 * @param img Image to be saved.
 * @param remain Pointer to the remaining list.
 * @param tail Pointer to the tail of the remaining list.
 */
void save_img_to_remaining_list(IMAGE_t img, IMG_DATA_t **remain, IMG_DATA_t **tail) {
    IMAGE_t *imgp = (IMAGE_t *) malloc(sizeof(IMAGE_t));
    *imgp = img;
    IMG_DATA_t *tmp = (IMG_DATA_t *) malloc(sizeof(IMG_DATA_t));
    tmp->to_image = imgp;
    tmp->next = NULL;
    if (!(*remain)) {
        *remain = tmp;
    } else {
        (*tail)->next = tmp;
    }
    *tail = tmp;
}

/**
 * @brief Generate element files
 *
 * Macro formally called "DUMBBELL" in the original code.
 *
 * loop over all_iprot
 *   Read in MSPS from file msp_file
 *   if even index than use mprot[iprot[i]->index/2]->pe as partner
 *   else use mprot[iprot[i]->index/2]->se as partner
 *   if ( iprot[i]->to_msp->score )
 *     if ( !ep_array[iprot[i]->ele_index-1]->flag )
 *       open ele_def_res/e# file
 *       print header
 *       set flag
 *      if ( partner_index != iprot[i]->ele_index )
 *        if ( iprot[i]->index%2 )
 *           print msp line [ partner, ele_index ] order
 *        else
 *           print msp line [ ele_index, partner ] order
 *      else
 *        if ( iprot[i]->index%2 )
 *          print msp line [ partner, ele_index ] order
 *
 * @param iprot_ct Count of IPROT_t structures.
 * @param all_iprot Array of IPROT_t structures.
 * @param msp_file File pointer to MSP file.
 * @param all_mprot Array of MSP prototypes.
 * @param ep_array Array of element prototypes.
 * @param all_ele Pointer to output file pointer for elements (to be updated).
 * @param err Error file pointer.
 * @param size_list Output file pointer for size list.
 * @param seq_names Array of sequence names.
 * @param seq_no Number of sequences.
 */
void build_element_files(int iprot_ct, IPROT_t **all_iprot, IPROT_t **iprot_shadow, FILE *msp_file, MPROT_t **all_mprot, EPROT_t **ep_array, FILE **all_ele, FILE *err, FILE *size_list, char **seq_names, int seq_no) {
    int i, partner_index;
    char ele_name[50];

    img_charge(iprot_shadow, iprot_ct, msp_file, seq_names, seq_no);
    for (i = 0; i < iprot_ct; i++) {
        if ((*(all_iprot + i))->index % 2) {
            partner_index = (*(all_mprot + (*(all_iprot + i))->index / 2))->pe;
        } else {
            partner_index = (*(all_mprot + (*(all_iprot + i))->index / 2))->se;
        }
        if ((*(all_iprot + i))->to_msp->score) {
            if (!(*(ep_array + ((*(all_iprot + i))->ele_index - 1)))->flag) {
                sprintf(ele_name, "ele_def_res/e%d", (*(all_iprot + i))->ele_index);

                // Close the current all_ele file and open a new one
                if (*all_ele != NULL) {
                    fclose(*all_ele);
                }
                *all_ele = fopen(ele_name, "w");
                if (*all_ele == NULL) {
                    fprintf(err, "Error opening file: %s\n", ele_name);
                    continue; // Skip to next iteration if file can't be opened
                }

                fprintf(*all_ele, "index %d\n", (*(all_iprot + i))->ele_index);
                fprintf(*all_ele, "frag %s %ld %ld\n",
                        (*(ep_array + ((*(all_iprot + i))->ele_index - 1)))->frag.seq_name,
                        (*(ep_array + ((*(all_iprot + i))->ele_index - 1)))->frag.lb,
                        (*(ep_array + ((*(all_iprot + i))->ele_index - 1)))->frag.rb);
                fprintf(*all_ele, "img_no %d\n",
                        (*(ep_array + ((*(all_iprot + i))->ele_index - 1)))->img_no);
                fprintf(size_list, "%d %d\n",
                        (*(all_iprot + i))->ele_index,
                        (*(ep_array + ((*(all_iprot + i))->ele_index - 1)))->img_no);
                (*(ep_array + ((*(all_iprot + i))->ele_index - 1)))->flag = 1;
            }
            if (partner_index != (*(all_iprot + i))->ele_index) {
                if ((*(all_iprot + i))->index % 2) {
                    fprintf(*all_ele, "msp %d s %d %.1f %d %d %s %ld %ld %d %s %ld %ld\n",
                            (*(all_iprot + i))->index,
                            (*(all_iprot + i))->to_msp->score,
                            (*(all_iprot + i))->to_msp->iden,
                            (*(all_iprot + i))->to_msp->direction,
                            partner_index,
                            (*(all_iprot + i))->to_msp->query.frag.seq_name,
                            (*(all_iprot + i))->to_msp->query.frag.lb,
                            (*(all_iprot + i))->to_msp->query.frag.rb,
                            (*(all_iprot + i))->ele_index,
                            (*(all_iprot + i))->to_msp->sbjct.frag.seq_name,
                            (*(all_iprot + i))->to_msp->sbjct.frag.lb,
                            (*(all_iprot + i))->to_msp->sbjct.frag.rb);
                } else {
                    fprintf(*all_ele, "msp %d s %d %.1f %d %d %s %ld %ld %d %s %ld %ld\n",
                            (*(all_iprot + i))->index,
                            (*(all_iprot + i))->to_msp->score,
                            (*(all_iprot + i))->to_msp->iden,
                            (*(all_iprot + i))->to_msp->direction,
                            (*(all_iprot + i))->ele_index,
                            (*(all_iprot + i))->to_msp->query.frag.seq_name,
                            (*(all_iprot + i))->to_msp->query.frag.lb,
                            (*(all_iprot + i))->to_msp->query.frag.rb,
                            partner_index,
                            (*(all_iprot + i))->to_msp->sbjct.frag.seq_name,
                            (*(all_iprot + i))->to_msp->sbjct.frag.lb,
                            (*(all_iprot + i))->to_msp->sbjct.frag.rb);
                }
            } else {
                if ((*(all_iprot + i))->index % 2) {
                    fprintf(*all_ele, "msp %d s %d %.1f %d %d %s %ld %ld %d %s %ld %ld\n",
                            (*(all_iprot + i))->index,
                            (*(all_iprot + i))->to_msp->score,
                            (*(all_iprot + i))->to_msp->iden,
                            (*(all_iprot + i))->to_msp->direction,
                            partner_index,
                            (*(all_iprot + i))->to_msp->query.frag.seq_name,
                            (*(all_iprot + i))->to_msp->query.frag.lb,
                            (*(all_iprot + i))->to_msp->query.frag.rb,
                            (*(all_iprot + i))->ele_index,
                            (*(all_iprot + i))->to_msp->sbjct.frag.seq_name,
                            (*(all_iprot + i))->to_msp->sbjct.frag.lb,
                            (*(all_iprot + i))->to_msp->sbjct.frag.rb);
                }
            }
        } // If score
    } // For loop
}

MU_TEST(test_get_seq_index) { 
    // Setup test data
    int seq_no = 3;
    char ** seq_names = malloc(seq_no * sizeof(char *));
    seq_names[0] = "seq1";
    seq_names[1] = "seq2";
    seq_names[2] = "seq3";

    mu_assert(get_seq_index(0, 2, "seq2", seq_names) == 1, "error, seq_index != 1");

    free(seq_names);
}

// Additional unit test function examples
MU_TEST(test_sing_cov) {
    FRAG_t f1 = {"seq1", 0, 100};
    FRAG_t f2 = {"seq1", 50, 150};

    mu_assert(sing_cov(&f1, &f2, 0.5) == 1, "error, sing_cov != 1");
    mu_assert(sing_cov(&f1, &f2, 0.8) == 0, "error, sing_cov != 0");
}

MU_TEST(test_doub_cov) {
    FRAG_t f1 = {"seq1", 0, 100};
    FRAG_t f2 = {"seq1", 50, 150};

    mu_assert(doub_cov(&f1, &f2, 0.1) == 1, "error, doub_cov != 1");
    mu_assert(doub_cov(&f1, &f2, 0.8) == 0, "error, doub_cov != 0");
}

/**
 * @brief Unit test for get_seq_names function.
 */
MU_TEST(test_get_seq_names) {
    FILE *seq_list;
    char **seq_names = NULL;
    int seq_count;

    // Create a temporary file for testing
    seq_list = fopen("temp_seq_list.txt", "w+");
    mu_assert(seq_list, "error, could not create temp file for test_get_seq_names");

    // Write test data to the temporary file
    fprintf(seq_list, "3\nseq1\n>seq2\n> seq3\n");
    rewind(seq_list);

    // Run get_seq_names
    seq_count = get_seq_names(seq_list, &seq_names);

    mu_assert(seq_count == 3, "error, get_seq_names did not read correct count");
    mu_assert(strcmp(seq_names[0], "seq1") == 0, "error, get_seq_names failed to read first sequence");
    mu_assert(strcmp(seq_names[1], "seq2") == 0, "error, get_seq_names failed to read second sequence");
    mu_assert(strcmp(seq_names[2], "seq3") == 0, "error, get_seq_names failed to read third sequence");

    // Clean up
    fclose(seq_list);
    remove("temp_seq_list.txt");

    for (int i = 0; i < seq_count; i++) {
        free(seq_names[i]);
    }
    free(seq_names);
}

MU_TEST_SUITE(all_tests) {
    MU_RUN_TEST(test_get_seq_index);
    MU_RUN_TEST(test_sing_cov);
    MU_RUN_TEST(test_doub_cov);
    MU_RUN_TEST(test_get_seq_names);
}

void print_usage() {
    printf("Usage: eledef -s <seq_list> -m <msp_file> -t <method> [-c <cutoff>] [-o <min_overlap>] [--test]\n");
    printf("  -s, --seq_list <seq_list>       Input file for sequence names list\n");
    printf("  -m, --msp_file <msp_file>       Input file for MSPs\n");
    printf("  -t, --method <method>           Method for evaluating overlap ('single' or 'double')\n");
    printf("  -c, --cutoff <cutoff>           Optional: cutoff value (default 0.5 for single, 0.9 for double)\n");
    printf("  -o, --min_overlap <min_overlap> Optional: minimum overlap between fragments (default 10)\n");
    printf("  --test                          Run inline unit tests\n");
}

int main(int argc, char *argv[]) {
    int ele_ct = 0, msp_ct, i, method = 0;
    float cutoff = 0.0;
    char line[150];
    char *seq_list_filename = NULL;
    char *msp_file_filename = NULL;
    char *method_str = NULL;
    MPROT_t **all_mprot;
    EPROT_t *all_ep = NULL, **ep_array, *ep_tmp;
    IPROT_t **all_iprot, **iprot_shadow;
    int iprot_ct = 0, partner_index;
    FILE *frags, *seq_list, *msp_file;
    FILE *msp_no, *all_ele, *img_prot, *ele_no, *err, *size_list;
    int seq_no;
    char **seq_names;

    static struct option long_options[] = {
        {"seq_list", required_argument, 0, 's'},
        {"msp_file", required_argument, 0, 'm'},
        {"method", required_argument, 0, 't'},
        {"cutoff", optional_argument, 0, 'c'},
        {"test", no_argument, 0, 0},
        {0, 0, 0, 0}
    };

    int opt;
    while ((opt = getopt_long(argc, argv, "s:m:t:c:", long_options, NULL)) != -1) {
        switch (opt) {
            case 's':
                seq_list_filename = optarg;
                break;
            case 'm':
                msp_file_filename = optarg;
                break;
            case 't':
                method_str = optarg;
                if (strcmp(method_str, "single") == 0) {
                    method = 1;
                    cutoff = 0.5;
                } else if (strcmp(method_str, "double") == 0) {
                    method = 2;
                    cutoff = 0.9;
                } else {
                    fprintf(stderr, "Invalid method: %s. Choose 'single' or 'double'.\n", method_str);
                    print_usage();
                    exit(EXIT_FAILURE);
                }
                break;
            case 'c':
                cutoff = atof(optarg);
                break;
            case 0:
                // Run unit tests if '--test' is specified
                {
                    MU_RUN_SUITE(all_tests);
                    MU_REPORT();
                    return MU_EXIT_CODE;
                }
                break;
            default:
                print_usage();
                exit(EXIT_FAILURE);
        }
    }

    if (!seq_list_filename || !msp_file_filename || method == 0) {
        print_usage();
        exit(EXIT_FAILURE);
    }

    seq_list = fopen(seq_list_filename, "r");
    if (!seq_list) {
        fprintf(stderr, "Input file for sequence name list %s not found. Exit.\n", seq_list_filename);
        return 2;
    }

    seq_no = get_seq_names(seq_list, &seq_names);
    if ( seq_no <= 0 ) {
        fprintf(stderr, "Count of sequence names in %s < 1. Exit.\n", seq_list_filename);
        return 2;
    }

    msp_file = fopen(msp_file_filename, "r");
    if (!msp_file) {
        fprintf(stderr, "Input file of MSPs %s not found. Exit.\n", msp_file_filename);
        return 2;
    }

    /* Open other required files */
    if (!(frags = fopen("images/images_sorted", "r"))) {
        fprintf(stderr, "Cannot open the fragment list file, exiting.\n");
        exit(1);
    }
    if (!(msp_no = fopen("summary/ori_msp_no", "r"))) {
        fprintf(stderr, "Cannot open msp_no, exiting.\n");
        exit(1);
    }
    err = fopen("ele_def_res/errors", "w");
    all_ele = fopen("summary/naive_eles", "w");
    img_prot = fopen("ele_def_res/img_prot", "w");
    ele_no = fopen("summary/naive_ele_no", "w");
    size_list = fopen("ele_def_res/size_list", "w");

    while (fgets(line, 100, msp_no)) {
        msp_ct = atoi(line);
    }

    all_mprot = (MPROT_t **) malloc(msp_ct * sizeof(MPROT_t *));
    for (i = 0; i < msp_ct; i++) {
        all_mprot[i] = (MPROT_t *) malloc(sizeof(MPROT_t));
        all_mprot[i]->pe = 0;
        all_mprot[i]->se = 0;
    }

    // Initial clustering into elements
    ele_def(method, frags, img_prot, cutoff, &all_ep, &ele_ct, all_mprot, all_ele, ele_no, seq_names, seq_no);

    fclose(ele_no);
    fclose(img_prot);
    fclose(msp_no);
    fclose(frags);

    img_prot = fopen("ele_def_res/img_prot", "r");

    // Initialize structures to hold elements
    ep_array = (EPROT_t **) malloc(ele_ct * sizeof(EPROT_t *));
    ep_tmp = all_ep;
    i = 0;
    while (ep_tmp) {
        ep_array[i++] = ep_tmp;
        ep_tmp = ep_tmp->next;
    }

    // Initialize structures to hold images
    all_iprot = (IPROT_t **) malloc(IMG_CAP * sizeof(IPROT_t *));
    iprot_shadow = (IPROT_t **) malloc(IMG_CAP * sizeof(IPROT_t *));
    for (i = 0; i < IMG_CAP; i++) {
        all_iprot[i] = (IPROT_t *) malloc(sizeof(IPROT_t));
        all_iprot[i]->to_msp = (MSP_t *) malloc(sizeof(MSP_t));
        iprot_shadow[i] = all_iprot[i];
    }

    printf("Reading images...\n");
    // Read in Element/Image tuples and write out log files for elements
    while (fgets(line, 100, img_prot)) {
        sscanf(line, "%d %d\n", &all_iprot[iprot_ct]->ele_index, &all_iprot[iprot_ct]->index);
        iprot_ct++;
        if (iprot_ct == IMG_CAP) {
            // build element files for a batch of IMG_CAP images
            build_element_files(iprot_ct, all_iprot, iprot_shadow, msp_file, all_mprot, ep_array, &all_ele, err, size_list, seq_names, seq_no);
            iprot_ct = 0;
        }
    }

    // trailing case, if there was an incomplete batch at the end
    if (iprot_ct) {
        for (i = 0; i < iprot_ct; i++) {
            iprot_shadow[i] = all_iprot[i];
        }
        build_element_files(iprot_ct, all_iprot, iprot_shadow, msp_file, all_mprot, ep_array, &all_ele, err, size_list, seq_names, seq_no);
    }

    // Clean up memory
    for (i = 0; i < IMG_CAP; i++) {
        free(all_iprot[i]->to_msp);
        free(all_iprot[i]);
    }
    free(all_iprot);
    free(iprot_shadow);

    for (i = 0; i < msp_ct; i++) {
        free(all_mprot[i]);
    }
    free(all_mprot);

    for (i = 0; i < ele_ct; i++) {
        free(ep_array[i]);
    }
    free(ep_array);

    fclose(all_ele);
    fclose(err);
    fclose(size_list);

    return 0;
}
