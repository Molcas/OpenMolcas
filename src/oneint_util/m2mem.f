************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 1993, Roland Lindh                                     *
************************************************************************
      Subroutine M2Mem(nHer,MemM2,la,lb,lr)
************************************************************************
*  Object: to compute the number of real*8 the kernal routine will     *
*          need for the computation of a matrix element between two    *
*          cartesian Gaussin functions with the total angular momentum *
*          of la and lb (la=0 s-function, la=1 p-function, etc.)       *
*          lr is the order of the operator (this is only used when the *
*          integrals are computed with the Hermite-Gauss quadrature).  *
*                                                                      *
*  Called from: OneEl                                                  *
*                                                                      *
************************************************************************
*
      nElem(i) = (i+1)*(i+2)/2
*
      nHer=(la+lb+2)/2
      MemM2 = 3*nHer*(la+1) +
     &        3*nHer*(lb+1) +
     &        3*nHer +
     &        3*(la+1)*(lb+1) +
     &        5 + nElem(la)*nElem(lb)
*
      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_integer(lr)
      End
*
      Subroutine PAM2Mem(nHer,MemM2,la,lb,lr)
************************************************************************
*  Object: to compute the number of real*8 the kernal routine will     *
*          need for the computation of a matrix element between two    *
*          cartesian Gaussin functions with the total angular momentum *
*          of la and lb (la=0 s-function, la=1 p-function, etc.)       *
*          lr is the order of the operator (this is only used when the *
*          integrals are computed with the Hermite-Gauss quadrature).  *
*                                                                      *
*  Called from: OneEl                                                  *
*                                                                      *
************************************************************************
*
      nElem(i) = (i+1)*(i+2)/2
*
      nComp = nElem(lr)
*
       nHer=(la+lb+lr+2)/2
      MemM2 = 3*nHer*(la+1) +
     &        3*nHer*(lb+1) +
     &        3*nHer*(lr+1) +
     &        3*(la+1)*(lb+1)*(lr+1) +
     &        5 + nElem(la)*nElem(lb)*nComp
*
      Return
      End
