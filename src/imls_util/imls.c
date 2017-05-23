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


#include <stdio.h>
#include <stdlib.h>

#include <math.h>
#include <assert.h>
#include <string.h>

#include <imls.h>
#include <imls_plan.h>
#include <svd.h>

void init_imls(struct imls *fit, int dof, int m,
                       int nr, double **data, double *energy)
{
        int i, j, lbasis;

        fit->multi = 0;

        fit->plan = init_imls_plan(dof, m);
        print_plan(fit->plan, m);

        lbasis = fit->lbasis = get_plan_size(fit->plan, m);
        fit->nr = nr;
        fit->m = m;
        fit->dof = dof;

        fit->data = data;
        fit->energy = energy;

        /* allocate storage */

        /* weighting */
        fit->W    = (double*)malloc(sizeof(double)*nr);

        /* least-squares polynomial coefficients  */
        fit->a    = (double*)malloc(sizeof(double)*lbasis);
        fit->b    = (double*)malloc(sizeof(double)*lbasis);

        /* diagonal values from SVD  */
        fit->wsvd = (double*)malloc(sizeof(double)*lbasis);

        /* mononomials at evaluation point  */
        fit->basis = (double*)malloc(sizeof(double)*lbasis);

        /* for matrices */

        fit->B = (double**)malloc(sizeof(double*)*nr);
        fit->B[0] = (double*)malloc(sizeof(double)*nr*lbasis);

        fit->u = (double**)malloc(sizeof(double*)*lbasis);
        fit->u[0] = (double*)malloc(sizeof(double)*lbasis*lbasis);

        fit->v = (double**)malloc(sizeof(double*)*lbasis);
        fit->v[0] = (double*)malloc(sizeof(double)*lbasis*lbasis);

        fit->tmp = (double**)malloc(sizeof(double*)*lbasis);
        fit->tmp[0] = (double*)malloc(sizeof(double)*nr*lbasis);

        for (i=1; i<nr; i++) {
                fit->B[i] = fit->B[i-1] + lbasis;
        }

        for (j=1; j<lbasis; j++) {
                fit->u[j] = fit->u[j-1] + lbasis;
                fit->v[j] = fit->v[j-1] + lbasis;
                fit->tmp[j] = fit->tmp[j-1] + nr;
        }

        /* initialise basis*/
        for (j=0; j < lbasis; j++)
                fit->basis[j] = 0.0;

        /* fill B matrix with plan */
        for (i=0; i < nr; i++)
                apply_plan(fit->B[i], data[i], lbasis, fit->plan);

        #ifdef DEBUG
        printf("total number of elements in plan: %d \n",
                        get_plan_size(fit->plan,m));

        for (i=0; i < nr; i++)
                for (j=0; j < lbasis; j++)
                printf("B[%d][%d] = %10.8f \n", i, j, fit->B[i][j]);
        #endif

}

void init_imls_multi(struct imls *fit, double **grad)
{
        int i, j, k, dof, nr, m, lbasis;

        dof = fit->dof;
        m   = fit->m;
        nr  = fit->nr;
        lbasis = fit->lbasis;

        fit->multi = 1;

        /* alloc surface for every dof*/
        fit->surf = (struct multi*)malloc(sizeof(struct multi)*dof);

        for (k=0; k < dof; k++) {
                fit->surf[k].plan = init_imls_multi_plan(k, m, fit->plan);
                print_plan(fit->surf[k].plan, m);

                fit->surf[k].df = grad[k];

                fit->surf[k].db = (double*)malloc(sizeof(double)*lbasis);

                fit->surf[k].dB = (double**)malloc(sizeof(double*)*nr);
                fit->surf[k].dB[0] = (double*)malloc(sizeof(double)*nr*lbasis);

                fit->surf[k].du = (double**)malloc(sizeof(double*)*lbasis);
                fit->surf[k].du[0] = (double*)malloc(sizeof(double)*lbasis*lbasis);

                fit->surf[k].tmp = (double**)malloc(sizeof(double*)*lbasis);
                fit->surf[k].tmp[0] = (double*)malloc(sizeof(double)*nr*lbasis);

                for (i=1; i<nr; i++)
                        fit->surf[k].dB[i] = fit->surf[k].dB[i-1] + lbasis;

                for (j=1; j<lbasis; j++) {
                        fit->surf[k].du[j] = fit->surf[k].du[j-1] + lbasis;
                        fit->surf[k].tmp[j] = fit->surf[k].tmp[j-1] + nr;
                }

                /* fill B matrix with plan */
                for (i=0; i < nr; i++)
                        apply_plan(fit->surf[k].dB[i], fit->data[i], lbasis,
                                        fit->surf[k].plan);
                #ifdef DEBUG
                printf("total number of elements in plan: %d \n",
                                get_plan_size(fit->surf[k].plan, m));

                for (i=0; i < nr; i++)
                        for (j=0; j < lbasis; j++)
                        printf("B[%d][%d] = %10.8f \n", i, j, fit->surf[k].dB[i][j]);
                #endif
        }
}

void free_imls(struct imls *fit)
{
        int i, k;

        if (fit == NULL) {
                return;
        }

        if (fit->plan != NULL) {
                for (i = 0; i < (fit->m)+1; i++)
                        free_imls_plan(fit->plan[i]);
                free(fit->plan);
        }

        free(fit->wsvd);
        free(fit->W);
        free(fit->a);
        free(fit->b);
        free(fit->basis);
        free(fit->B[0]);
        free(fit->B);
        free(fit->u[0]);
        free(fit->u);
        free(fit->v[0]);
        free(fit->v);
        free(fit->tmp[0]);
        free(fit->tmp);

        if (fit->multi) {
                for (k=0; k < fit->dof; k++)
                        free_multi(&(fit->surf[k]), fit->m);
                free(fit->surf);
        }

}

void free_multi(struct multi *surf, int m)
{
        int i;
        if (surf->plan != NULL) {
                for (i = 0; i < (m+1); i++)
                        free_imls_plan(surf->plan[i]);
                free(surf->plan);
        }

        free(surf->db);
        free(surf->dB[0]);
        free(surf->dB);
        free(surf->du[0]);
        free(surf->du);
        free(surf->tmp[0]);
        free(surf->tmp);

}

double weighting(double *eval, double *data, int dof)
{
        int i;
        double length, length2 = 0.0, tmp = 0.0;

        for (i=0; i < dof; i++) {
                tmp = eval[i] - data[i];
                length2 += tmp*tmp;
        }

        length = sqrt(length2);

        return exp(-(length*length))/(pow(length,6)+1.0e-10);
}

void imls_execute(struct imls *fit, double *eval_point)
{
        int i, j, k;
        double interpolated = 0.0;

        if (fit->multi) {
                imls_multi_execute(fit, eval_point);
                return;
        }

        /* populate weighting vector */
        for (i=0; i < fit->nr; i++)
                fit->W[i] = weighting(eval_point, fit->data[i], fit->dof);

        /* tmp = B^t * W */
        for (j=0; j < fit->lbasis; j++)
                for (i=0; i < fit->nr; i++)
                        fit->tmp[j][i] = fit->B[i][j] * fit->W[i];

        /* u = B^t * W * B = tmp * B */
        memset(fit->u[0], '\0', sizeof(fit->u[0][0])*(fit->lbasis)*(fit->lbasis));

        for (k=0; k < fit->lbasis; k++)
                for (j=0; j < fit->lbasis; j++)
                        for (i=0; i < fit->nr; i++)
                                fit->u[k][j] += fit->tmp[k][i] * fit->B[i][j];

        /* Singular Value Decomposition */
        svdcmp(fit->u, fit->lbasis, fit->lbasis, fit->wsvd, fit->v);

        /* b = B^t * W  * f = tmp * f */
        for (j=0; j < fit->lbasis; j++)
                fit->b[j] = 0.0;

        for (i=0; i < fit->nr; i++)
                for (j=0; j < fit->lbasis; j++)
                        fit->b[j] += fit->tmp[j][i] * fit->energy[i];

        /* Calculate coefficients for polynomial */
        svbksb(fit->u, fit->lbasis, fit->lbasis, fit->wsvd,
                           fit->v, fit->b, fit->a);

        /* interpolated values */
        apply_plan(fit->basis, eval_point, fit->lbasis, fit->plan);

        for (j=0; j < fit->lbasis; j++)
                interpolated += fit->a[j] * fit->basis[j];

        #ifdef DEBUG
        for (j=0; j < fit->lbasis; j++)
                printf("basis[%d] = %10.8f \n", j, fit->basis[j]);
        for (j=0; j < fit->lbasis; j++)
                printf("a[%d] = %10.8f \n", j, fit->a[j]);


        printf("INTERPOLATED VALUE = %10.6f \n", interpolated);
        #endif

}

/* MULTI-SURFACE EXECUTION */

void imls_multi_execute(struct imls *fit, double *eval_point)
{
        int i, j, t, k;
        double interpolated = 0.0;

        /* populate weighting vector */
        for (i=0; i < fit->nr; i++)
                fit->W[i] = weighting(eval_point, fit->data[i], fit->dof);

        /* tmp = B^t * W */
        for (j=0; j < fit->lbasis; j++)
                for (i=0; i < fit->nr; i++)
                        fit->tmp[j][i] = fit->B[i][j] * fit->W[i];

        /* and for multi-surfaces */
        for (k=0; k < fit->dof; k++)
                for (j=0; j < fit->lbasis; j++)
                        for (i=0; i < fit->nr; i++)
                                fit->surf[k].tmp[j][i] =
                                        fit->surf[k].dB[i][j] * fit->W[i];

        /* u = B^t * W * B = tmp * B */
        for (t=0; t < fit->lbasis; t++)
                for (j=0; j < fit->lbasis; j++)
                        fit->u[t][j] = 0.0;

        for (t=0; t < fit->lbasis; t++)
                for (j=0; j < fit->lbasis; j++)
                        for (i=0; i < fit->nr; i++)
                                fit->u[t][j] += fit->tmp[t][i] * fit->B[i][j];

        /* and for multi-surfaces */
        for (k=0; k < fit->dof; k++)
                for (t=0; t < fit->lbasis; t++)
                        for (j=0; j < fit->lbasis; j++)
                                fit->surf[k].du[k][j] = 0.0;

        for (k=0; k < fit->dof; k++)
                for (t=0; t < fit->lbasis; t++)
                        for (j=0; j < fit->lbasis; j++)
                                for (i=0; i < fit->nr; i++)
                                        fit->surf[k].du[t][j] +=
                                                fit->surf[k].tmp[t][i] * fit->surf[k].dB[i][j];

        /* b = B^t * W  * f = tmp * f */
        for (j=0; j < fit->lbasis; j++)
                fit->b[j] = 0.0;

        for (i=0; i < fit->nr; i++)
                for (j=0; j < fit->lbasis; j++)
                        fit->b[j] += fit->tmp[j][i] * fit->energy[i];

        /* and for multi-surfaces */
        for (k=0; k < fit->dof; k++)
                for (j=0; j < fit->lbasis; j++)
                        fit->surf[k].db[j] = 0.0;

        for (k=0; k < fit->dof; k++)
                for (i=0; i < fit->nr; i++)
                        for (j=0; j < fit->lbasis; j++)
                                fit->surf[k].db[j] +=
                                        fit->surf[k].tmp[j][i] * fit->surf[k].df[i];

        /* calculate lambda :
           relative control variable for multi surface fitting */
        fit->lambda = 1.0;
        for (k=0; k < fit->dof; k++)
                fit->surf[k].lambda = 1.0;

        /* merge multi-surfaces together: u & b */
        for (t=0; t < fit->lbasis; t++)
                for (j=0; j < fit->lbasis; j++) {
                        fit->u[t][j] *= fit->lambda * fit->lambda;
                        for (k=0; k < fit->dof; k++)
                                fit->u[t][j] += fit->surf[k].lambda * fit->surf[k].lambda *
                                        fit->surf[k].du[t][j];
                }

        for (j=0; j < fit->lbasis; j++) {
                fit->b[j] *= fit->lambda * fit->lambda;
                for (k=0; k < fit->dof; k++)
                        fit->b[j] += fit->surf[k].lambda * fit->surf[k].lambda *
                                                        fit->surf[k].db[j];
        }

        /* Singular Value Decomposition */
        svdcmp(fit->u, fit->lbasis, fit->lbasis, fit->wsvd, fit->v);

        /* Calculate coefficients for polynomial */
        svbksb(fit->u, fit->lbasis, fit->lbasis, fit->wsvd,
                           fit->v, fit->b, fit->a);

        /* interpolated values */
        apply_plan(fit->basis, eval_point, fit->lbasis, fit->plan);

        for (j=0; j < fit->lbasis; j++)
                interpolated += fit->a[j] * fit->basis[j];

        #ifdef DEBUG
        for (j=0; j < fit->lbasis; j++)
                printf("basis[%d] = %10.8f \n", j, fit->basis[j]);
        for (j=0; j < fit->lbasis; j++)
                printf("a[%d] = %10.8f \n", j, fit->a[j]);


        printf("INTERPOLATED VALUE = %10.6f \n", interpolated);

        for (k=0; k < fit->dof; k++) {
                apply_plan(fit->basis, eval_point, fit->lbasis, fit->surf[k].plan);

                interpolated = 0.0;
                for (j=0; j < fit->lbasis; j++)
                        interpolated += fit->a[j] * fit->basis[j];
                printf("GRADIENT %d VALUE = %10.6f \n", k, interpolated);
        for (j=0; j < fit->lbasis; j++)
                printf("gradient basis[%d] = %10.8f \n", j, fit->basis[j]);
        }
        #endif

}

/* Relax on imls surface */

void imls_relax(struct imls *fit, double *eval_point)
{
        int i, k, dof, lbasis, iter;
        int max_iter, max_int, int_count;

        double tmp, energy, energy_last;
        double abs2, abs2_last, threshold, scale;
        double max_scale, min_scale, scl_ratio;

        double *grad_current = NULL;
        double *eval_last = NULL;

        /* initialization */
        dof              = fit->dof;
        lbasis    = fit->lbasis;
        threshold = 1.0e-8;
        max_scale = 1.0;
        min_scale = 0.01;
        scl_ratio = 8.0/10.0;
        scale     = 0.4;
        iter          = 1;
        max_iter  = 200;
        max_int   = 100;
        abs2_last   = abs2   = 0.0;
        energy_last = energy = 0.0;

        /* allocate memory */
        grad_current = (double*)malloc(sizeof(double)*dof);
        eval_last     = (double*)malloc(sizeof(double)*dof);

        for (i=1; i < dof; i++) {
                grad_current[i] = 0.0;
                eval_last[i] = eval_point[i];
        }

        /* realx */
        do {

                /* determine scaling */
                if (iter > 1) {
                        if (abs2_last < abs2) {
                                scale *= scl_ratio;
                        } else {
                                scale /= scl_ratio;
                        }
                        if (scale < min_scale)
                                scale = min_scale;
                        else if (scale > max_scale)
                                scale = max_scale;
                }

                /* move along gradient */
                for (i=0; i < dof; i++) {
                        eval_point[i] += scale * grad_current[i];
                }

                tmp = scale;
                int_count = 0;
                while(int_count < max_int) {

                /* perform new evaluation */
                imls_execute(fit, eval_point);

                /* calculate energy */
                apply_plan(fit->basis, eval_point, fit->lbasis, fit->plan);

                energy = 0.0;
                for (i=0; i < fit->lbasis; i++)
                        energy += fit->a[i] * fit->basis[i];

                /* reject step if energy increases */
                if ( iter > 1 && energy > energy_last) {
                        tmp *= 0.9;
                        for (i=0; i < dof; i++) {
                                eval_point[i] = eval_last[i] + tmp * grad_current[i];
                        }
                        int_count++;
                        continue;
                }
                break;
                }

                /* calculate new gradient */
                abs2 = 0.0;
                for (i=0; i < dof; i++) {
                        /* recipe for derivative */
                        apply_plan(fit->basis, eval_point, fit->lbasis,
                                           fit->surf[i].plan);

                        tmp = 0.0;
                        for (k=0; k < lbasis; k++)
                                tmp += fit->a[k] * fit->basis[k];

                        abs2 += tmp * tmp;
                        grad_current[i] = -tmp;
                }

                #ifdef DEBUG
                printf("****************************\n");
                printf("*** R E L A X   S T A T ****\n");
                printf("****************************\n");
                printf("Iteration           %8d\n", iter);
                printf("Thresh          %8.4f\n", threshold);
                printf("Abs**2          %8.4f\n", abs2);
                printf("Scale           %8.4f\n", scale);
                printf("Gradient\n");
                for (i=0; i < dof; i++)
                        printf("%6d          %8.4f\n", i, grad_current[i]);
                printf("Structure\n");
                for (i=0; i < dof; i++)
                        printf("%6d          %8.4f\n", i, eval_point[i]);
                #endif

                if (iter > max_iter)
                        /* not converged */
                        break;

                for (i=0; i < dof; i++)
                        eval_last[i] = eval_point[i];

                energy_last = energy;
                abs2_last = abs2;

                iter++;

        } while(abs2 > threshold );

        /* free memory */
        free(grad_current);

}
