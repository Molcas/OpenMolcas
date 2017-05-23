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
      Subroutine M1MmG(nRys,MmM1g,la,lb,lr)
************************************************************************
*                                                                      *
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
      Integer iAng(4)
*
      nElem(i) = (i+1)*(i+2)/2
*
      iAng(1) = la
      iAng(2) = lb
      iAng(3) = 0
      iAng(4) = 0
      Call MemRg1(iAng,nRys,Mem)
      MmM1g= 8 + Mem + nElem(la)*nElem(lb)
*
      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_integer(lr)
      End
