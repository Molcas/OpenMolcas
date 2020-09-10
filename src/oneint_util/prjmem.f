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
      Subroutine PrjMem(nHer,MemPrj,la,lb,lr)
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
      use Basis_Info, only: dbsc, nCnttp, Shells
#include "itmax.fh"
#include "info.fh"
*
      nElem(i) = (i+1)*(i+2)/2
*
      nHer=0
      MemPrj = 0
      Do 1960 iCnttp = 1, nCnttp
         If (.Not.dbsc(iCnttp)%ECP) Cycle
         Do 1966 iAng = 0, dbsc(iCnttp)%nPrj-1
            iShll = dbsc(iCnttp)%iPrj + iAng
            nExpi=Shells(iShll)%nExp
            nBasisi=Shells(iShll)%nBasis
            If (nExpi.eq.0 .or. nBasisi.eq.0) Go To 1966
*
            ip = 0
            nac = nElem(la)*nElem(iAng)
            ip = ip + nExpi*nac
            ip = ip + 3 * nExpi
            ip = ip + nExpi
            ip = ip + nExpi
            ip = ip + nExpi
*
            Call MltMmP(nH,MemMlt,la,iAng,lr)
            nHer = Max(nH,nHer)
            MemPrj = Max(MemPrj,ip+nExpi*MemMlt)
            ip = ip - 6 * nExpi
*
            ncb = nElem(iAng)*nElem(lb)
            ip = ip + nExpi*ncb
            ip = ip + 3 * nExpi
            ip = ip + nExpi
            ip = ip + nExpi
            ip = ip + nExpi
*
            Call MltMmP(nH,MemMlt,iAng,lb,lr)
            nHer = Max(nH,nHer)
            MemPrj = Max(MemPrj,ip+nExpi*MemMlt)
            ip = ip - 6 * nExpi
*
            ip = ip + Max(nExpi*nac,ncb*nBasisi)
            MemPrj = Max(MemPrj,ip)
*
 1966    Continue
 1960 Continue
*
      Return
      End
