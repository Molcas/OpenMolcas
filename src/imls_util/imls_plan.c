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
 *  ETH Zurich, Lund University
 *
 ***/

#ifdef DEBUG
#include <stdio.h>
#endif
#include <stdlib.h>
#include <assert.h>

#include <imls_plan.h>

struct plan_node* alloc_node()
        /*
         Allocating storage for a single node

         Returns: pointer to struct plan_node
         */
{
        struct plan_node *tmp =
                (struct plan_node *)malloc(sizeof(struct plan_node));

        tmp->val = NULL;
        tmp->dof = NULL;

        return(tmp);
}

void add_node(struct plan_node **start, int i, struct plan_node *node)
        /*
         Adding a node to plan structure

         In: start        pointer to plan structure
                 i                array indeces (= order of monomials)
                 node         node which is addded to plan
         */
{
        struct plan_node *tmp = NULL;

        if (*(start+i) == NULL)
                *(start+i) = node;
        else {

                tmp = *(start+i);

                while (tmp->next != NULL)
                        tmp = tmp->next;

                tmp->next = node;
        }
}

void free_imls_plan(struct plan_node *node)
        /*
         freeing nodes recursively

         In: node         node which freed
         */
{
        if (node == NULL)
                return;

        if (node->next != NULL)
                /* RECURSIVE DELETION */
                free_imls_plan(node->next);

        if (node->val != NULL)
                free(node->val);

        if (node->dof != NULL)
                free(node->dof);

        free(node);

}

struct plan_node** init_imls_plan(int dof, int m)
        /*
         initialises plan for struct fit with
         a given number of DOF and a given order
         of polynomial m

         Returns: dof         number of degrees of freedom
                          m                order of desired polynomial
         */
{
        int i, j, k, kk;

        struct plan_node *above = NULL, *cur = NULL;

        /* node vector */
        struct plan_node **plan = NULL;
    plan = (struct plan_node **)malloc(sizeof(struct plan_node*)*(m+1));

        for (i=0; i<(m+1); i++)
                plan[i] = NULL;

        /* zeroth order */
        cur = *plan = alloc_node();
        cur->next = NULL;

        cur->val_length = 1;
        cur->val = (int*)malloc(sizeof(int)*(cur->val_length));
        cur->val[0] = -1;

        cur->dof_length = dof;
        cur->dof = (int*)malloc(sizeof(int)*(cur->dof_length));

        for (i=0; i < cur->dof_length; i++)
                cur->dof[i] = i;

        /* generating monomials up to desired order */
        for (i=1; i<(m+1); i++) {
                above = plan[i-1];

                do {
                        /* perform loop over dof of above */
                        assert(above != NULL);

                        int above_dof_length = above->dof_length;

                        for (j=0; j < above_dof_length; j++) {

                                cur = alloc_node();
                                cur->next = NULL;

                                /* update and add val */
                                cur->val_length = above->val_length + 1;
                                cur->val = (int*)malloc(sizeof(int)*(cur->val_length));
                                for (k=0; k < above->val_length; k++)
                                        cur->val[k] = above->val[k];
                                cur->val[above->val_length] = above->dof[j];

                                /* update dof */
                                cur->dof_length = above_dof_length - j;
                                cur->dof = (int*)malloc(sizeof(int)*(above_dof_length-j));
                                for (k=j, kk=0; k < above_dof_length; k++, kk++)
                                        cur->dof[kk] = above->dof[k];

                           /* add new node */
                           add_node(plan, i, cur);
                        }

                } while ((above = above->next) != NULL);
        }
        return plan;
}

struct plan_node** init_imls_multi_plan(int k, int m, struct plan_node **plan)
        /*
         calculating the derivative of polynomial
         stored in struct plan_node **plan.

         In: k                 degree of freedom (= dE/dk)
                 m                order of desired polynomial
                 plan        plan structure
         Returns: pointer to plan sturcuture
         */
{
        int i, j, t;
        struct plan_node *node = NULL, *ptr_to_plan = NULL;
        int count;

        /* node vector */
        struct plan_node **surf =
                (struct plan_node **)malloc(sizeof(struct plan_node*)*(m+1));

        for (i=0; i<(m+1); i++)
                surf[i] = NULL;

        for (i=0; i<(m+1); i++) {

                ptr_to_plan = plan[i];

                while (ptr_to_plan != NULL) {
                        node = alloc_node();
                        node->next = NULL;

                        /* calc grad to k */
                        node->val_length = ptr_to_plan->val_length;
                        node->val = (int*)malloc(sizeof(int)*(ptr_to_plan->val_length));

                        for (j=0, count=0, t=0; j < node->val_length; j++)
                                if (ptr_to_plan->val[j] == k)
                                        count++;
                                else
                                        node->val[t++] = ptr_to_plan->val[j];

                        if (count > 0) {
                                for (j=0; j < (count-1); j++)
                                                node->val[t++] = k;
                                node->val[t] = count*1000;
                        } else
                                for (j=0; j < node->val_length; j++)
                                        node->val[j] = -2;

                        add_node(surf, i, node);

                        ptr_to_plan = ptr_to_plan->next;
                }
        }

        return surf;
}


void apply_plan(double *res, double *data, int lbasis, struct plan_node **plan)
{
        int i, j, k;
        struct plan_node *node = NULL;
        double tmp;

        node = *plan;
        for (i=0, k=0; i < lbasis; i++) {

                for (j=0, tmp=1.0; j < node->val_length; j++) {
                        if (node->val[j] >=1000)
                                   tmp *= node->val[j]/1000;
                        else if (node->val[j] == -2)
                                   tmp = 0.0;
                        else if (node->val[j] == -1)
                                   tmp *= 1.0;
                        else
                                   tmp *= data[node->val[j]];
                }

                res[i] = tmp;

                if (node->next == NULL)
                        node = plan[++k];
                else
                        node = node->next;
        }
}

int get_plan_size(struct plan_node **plan, int m)
{
        int i;
        int ret = 0;
        struct plan_node *ptr_to_plan = NULL;

        if (plan == NULL)
                return 0;

        for (i=0; i< m+1; i++) {

                ptr_to_plan = plan[i];

                while (ptr_to_plan != NULL) {
                        ret++;
                        ptr_to_plan = ptr_to_plan->next;
                }

        }

        return(ret);

}

void print_node(struct plan_node *node)
{
        #ifdef DEBUG
        int j;

        printf("VALUES: ");

        for (j=0; j<node->val_length; j++)
                printf("%d   ", node->val[j]);

        printf("DOF: ");

        for (j=0; j<node->dof_length; j++)
                printf("%d   ", node->dof[j]);

        printf("\n");
        #endif
}

void print_plan(struct plan_node **plan, int m)
{
        #ifdef DEBUG
        int i;

        struct plan_node *node = NULL;

        if (plan == NULL) {
                printf("plan is NULL.\n");
                return;
        }

        printf("performing NULL check on plan\n");
        for (i=0; i<(m+1); i++) {
                printf("%d    :",i);
                if (plan[i] == NULL)
                        printf(" NULL \n");
                else
                        printf(" ok \n");

        }

        for (i=0; i<(m+1); i++) {

                if ((node = plan[i]) == NULL) {
                        continue;
                }

                printf("ORDER: %d \n", i);

                do {
                        print_node(node);
                } while ((node=node->next) != NULL);{

                }

        }
        #endif
}
