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
* Copyright (C) 2008,2009, Koni Marti                                  *
***********************************************************************/
/***
 *
 *  Koni Marti
 *  27.2.2008
 *
 ***/

#ifndef SVD_H
#define SVD_H

void svbksb(double **u, int m, int n, double *w,
            double **v, double *b, double *x);

void svdcmp(double **a, int m, int n,  double *w, double **v);

#endif
