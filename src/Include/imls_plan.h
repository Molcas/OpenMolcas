/***********************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 2009, Koni Marti                                       *
***********************************************************************/
/***
 *
 *  Koni Marti
 *
 *  04.02.2009
 *
 ***/

#ifndef IMLS_PLAN_H
#define IMLS_PLAN_H

/* plan nodes */
struct plan_node
{
        /* array storing monomials*/
        int *val;
        /* array length */
        int val_length;

        /* array storing remaining DOF */
        int *dof;
        /* array length */
        int dof_length;

        /* pointer to adjacent node */
        struct plan_node *next;
};

/* allocate memory for a single node */
struct plan_node* alloc_node();

/* initializes plan for polynomial */
struct plan_node** init_imls_plan(int dof, int m);

/* initializing derivative plan from plan */
struct plan_node** init_imls_multi_plan(int k, int m, struct plan_node **plan);

/* add node to existing plan */
void add_node(struct plan_node **start, int i, struct plan_node *node);

/* freeing allocated row of nodes recursively */
void free_imls_plan(struct plan_node *node);

/* applying plan to data and stores into res */
void apply_plan(double *res, double *data, int lbasis,
                        struct plan_node **plan);

/* traversing plan structure and counts nodes */
int         get_plan_size(struct plan_node **plan, int m);

/* print entire plan structure */
void print_plan(struct plan_node **plan, int m);

/* print content of single node */
void print_node(struct plan_node *node);

#endif
