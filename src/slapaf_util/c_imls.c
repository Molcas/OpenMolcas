/***********************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
***********************************************************************/
#include <stdio.h>
#include <stdlib.h>

#include <imls.h>
#include <molcastype.h>

#ifdef _CAPITALS_
#define imls_driver IMLS_DRIVER
#else
#ifndef ADD_
#define imls_driver imls_driver_
#endif
#endif

void imls_driver(INT *ptr_nr, INT *ptr_dof, double Energy[],
                     double g[], double q[], double dq[])

   /*
   A new geometry is determined by finding a stationary point on
   the potential energy surface constructed by the IMLS algorithm.

   from
   Haptic Quantum Chemistry,
   K. H. Marti and M. Reiher,
   J. Comput. Chem. (2009) DOI: 10.1002/jcc.21201
   */

/*      nr          : iteration counter                                 *
 *      dof         : total number of internal coordinates              *
 *      Energy      : the energy of each iteration                      *
 *      g(*,kIter)  : the gradient in the internal coordinates          *
 *      q(*,kIter)  : the internal coordinates                          *
 *      dq(*,kIter) : the internal coordinates                          *
 *                                                                      *
 *  To Access Gradient Information:                                     *
 *                                                                      *
 *  for (i=0; i < *kIter; i++)                                          *
 *    for (j=0; j < *nInter; j++)                                       *
 *      printf("%-12d%-12d%10.8f\n", i, j, g[i*(*nInter)+j]);           *
 *                                                                      *
 *                                                                      */


{
        int i, j, m;
        int nr, dof;

        double *eval_point = NULL;
        double **derivative = NULL;
        double **coord    = NULL;

        struct imls fit;

        dof = *ptr_dof;
        nr  = *ptr_nr;
        m   = (nr<4)?(nr-1):(3);
        if (m < 1) m=1;

        /* allocating memory for IMLS */
        eval_point = (double*)malloc(sizeof(double)*dof);

        derivative = (double**)malloc(sizeof(double *) * dof);
        derivative[0]= (double*)malloc(sizeof(double) * nr * dof);
        coord      = (double**)malloc(sizeof(double *) * nr);
        coord[0]   = (double*)malloc(sizeof(double) * nr * dof);

        for (i=1; i < dof; i++) {
                derivative[i] = derivative[i-1] + nr;
        }

        for (i=1; i < nr; i++) {
                coord[i]    = coord[i-1] + dof;
        }

        /* initialization of data */
        for (i=0; i < nr; i++) {
                for (j=0; j < dof; j++) {
                         coord[i][j]    = q[i*dof+j];
                        derivative[j][i] = -g[i*dof+j];
                }
        }
        for (i=0; i < dof; i++) {
                eval_point[i] = q[(nr-1)*dof+i];
        }

        printf("*******  K o n i   M a r t i  ********\n\n");

        /* initialize single-surface imls structure */
        init_imls(&fit, dof, m, nr, coord, Energy);

        /* initialize multi-surface structure */
        init_imls_multi(&fit, derivative);

        /* relax geometry on imls surface */
        imls_relax(&fit, eval_point);

        /* calculate shifts (dq) */
        for (i=0; i < dof; i++) {
                dq[(nr-1)*dof+i] =  eval_point[i] - q[(nr-1)*dof+i];
        }

        /* copy new geometry to q+1 position */
        for (i=0; i < dof; i++) {
                q[(nr-1)*dof+i] = eval_point[i];
        }

        /* freeing memory */
        free(eval_point);
        free(derivative[0]);
        free(derivative);
        free(coord[0]);
        free(coord);

        /* freeing imls structure and internal memory */
        free_imls(&fit);

}
