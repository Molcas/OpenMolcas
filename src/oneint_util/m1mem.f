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
      Subroutine M1Mem(nRys,MemM1,la,lb,lr)
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
      Integer iAng(4)
*
      nabSz(ixyz) = (ixyz+1)*(ixyz+2)*(ixyz+3)/6  - 1
*
      Call mHrr(la,lb,nFlop,nMem)
*
      iAng(1) = la
      iAng(2) = lb
      iAng(3) = 0
      iAng(4) = 0
      Call MemRys(iAng,MemPrm)
      MemM10= 6 + MemPrm
      nRys = (la+lb+2)/2
*
      k = nabSz(la+lb) - nabSz(Max(la,lb)-1)
*
*-----nMem : memory for Hrr
*     k    : memory for primitives in M1Int0
*     MemM10: scratch in M1Int0
*
      MemM1 = MemM10 + Max(nMem,k)
*
      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_integer(lr)
      End
