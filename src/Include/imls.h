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

#ifndef IMLS_H
#define IMLS_H

#include <imls_plan.h>

struct imls
{
        /* construction plan for polynom */
        struct plan_node **plan;

        /* multi-surface fitting */
        struct multi *surf;

        /* number of terms in plan */
        int lbasis;

        /* number of data points */
        int nr;

        /* order of polynomial */
        int m;

        /* degrees of freedom */
        int dof;

        /* flag for multi-surface (1 or 0) */
        int multi;

        /* lambda value for energy surface */
        double lambda;

        /* data points (N x dof) */
        double **data;

        /* energy values (N) */
        double *energy;

        /* SVD matrices */
        double **u, **v, *wsvd;

        /* storing basis for all points */
        double **B, **tmp;

        /* weight, coefficient, polynomial basis for evaluation */
        double *W, *a, *b, *basis;

};

struct multi
{
        /* construction plan for surface (gradient) */
        struct plan_node **plan;

        /* lambda value for gradient surface */
        double lambda;

        /* derivatives of total energy*/
        double *df;

        /* auxilliary matrices */
        double **dB, **du, **tmp, *db;

};

/* initializing struct imls,
   allocate memory and populate B matrix */
void init_imls(struct imls *fit, int dof, int m, int nr,
                           double **data, double *energy);

/* initializing struct multi,
   allocate memory and populate dB matrix */
void init_imls_multi(struct imls *fit, double **grad);

/* freeing allocate memory for struct imls */
void free_imls(struct imls *fit);

/* freeing allocate memory for struct multi */
void free_multi(struct multi *surf, int m);

/* performing IMLS fit for a single surface */
void imls_execute(struct imls *fit, double *eval_point);

/* performing IMLS fit for a multi surfaces */
void imls_multi_execute(struct imls *fit, double *eval_point);

/* relax geometry to stationary point  */
void imls_relax(struct imls *fit, double *eval_point);

/* calculate the weighting factor in error functional */
double weighting(double *eval, double *data, int dof);

#endif
