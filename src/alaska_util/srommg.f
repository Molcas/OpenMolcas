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
      Subroutine SROMmG(nHer,MmSROG,la,lb,lr)
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
#include "itmax.fh"
#include "info.fh"
*
      nElem(i) = (i+1)*(i+2)/2
*
      MmSROG = 0
      nOrder=0
      Do 1960 iCnttp = 1, nCnttp
         If (.Not.ECP(iCnttp)) Go To 1960
         Do 1966 iAng = 0, nSRO_Shells(iCnttp)-1
            iShll = ipSRO(iCnttp) + iAng
            If (nExp(iShll).eq.0) Go To 1966
*
            ip = 0
            ip = ip + nExp(iShll)**2
            nac = 4*nElem(la)*nElem(iAng)
            ip = ip + nExp(iShll)*nac
            ip = ip + 3 * nExp(iShll)
            ip = ip + nExp(iShll)
            ip = ip + nExp(iShll)
            ip = ip + nExp(iShll)
            nHer = ((la+1)+iAng+2)/2
            nOrder = Max(nHer,nOrder)
            ip = ip + nExp(iShll)*3*nHer*(la+2)
            ip = ip + nExp(iShll)*3*nHer*(iAng+1)
            ip = ip + nExp(iShll)*3*nHer*(lr+1)
            ip = ip + nExp(iShll)*3*nHer*(la+2)*(iAng+1)*(lr+1)
            ip = ip + nExp(iShll)
*
            MmSROG = Max(MmSROG,ip)
            ip = ip - nExp(iShll)
     &         * (6 + 3*nHer*((la+2) + (iAng+1) + (lr+1)
     &         +  (la+2)*(iAng+1)*(lr+1)) + 1)
*
            ncb = 4*nElem(iAng)*nElem(lb)
            ip = ip + nExp(iShll)*ncb
            ip = ip + 3 * nExp(iShll)
            ip = ip + nExp(iShll)
            ip = ip + nExp(iShll)
            ip = ip + nExp(iShll)
            nHer = ((lb+1)+iAng+2)/2
            nOrder = Max(nHer,nOrder)
            ip = ip + nExp(iShll)*3*nHer*(lb+2)
            ip = ip + nExp(iShll)*3*nHer*(iAng+1)
            ip = ip + nExp(iShll)*3*nHer*(lr+1)
            ip = ip + nExp(iShll)*3*nHer*(lb+2)*(iAng+1)*(lr+1)
            ip = ip + nExp(iShll)
*
            MmSROG = Max(MmSROG,ip)
            ip = ip - nExp(iShll)
     &         * (6 + 3*nHer*((lb+2) + (iAng+1) + (lr+1)
     &         +  (lb+2)*(iAng+1)*(lr+1)) + 1)
*
            ip = ip + Max(nExp(iShll)*nac,ncb*nExp(iShll))
            MmSROG = Max(MmSROG,ip)
*
 1966    Continue
 1960 Continue
*
      Return
      End
