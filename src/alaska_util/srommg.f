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
      use Basis_Info, only: nCnttp, Shells
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
            nExpi=Shells(iShll)%nExp
            If (nExpi.eq.0) Go To 1966
*
            ip = 0
            ip = ip + nExpi**2
            nac = 4*nElem(la)*nElem(iAng)
            ip = ip + nExpi*nac
            ip = ip + 3 * nExpi
            ip = ip + nExpi
            ip = ip + nExpi
            ip = ip + nExpi
            nHer = ((la+1)+iAng+2)/2
            nOrder = Max(nHer,nOrder)
            ip = ip + nExpi*3*nHer*(la+2)
            ip = ip + nExpi*3*nHer*(iAng+1)
            ip = ip + nExpi*3*nHer*(lr+1)
            ip = ip + nExpi*3*nHer*(la+2)*(iAng+1)*(lr+1)
            ip = ip + nExpi
*
            MmSROG = Max(MmSROG,ip)
            ip = ip - nExpi
     &         * (6 + 3*nHer*((la+2) + (iAng+1) + (lr+1)
     &         +  (la+2)*(iAng+1)*(lr+1)) + 1)
*
            ncb = 4*nElem(iAng)*nElem(lb)
            ip = ip + nExpi*ncb
            ip = ip + 3 * nExpi
            ip = ip + nExpi
            ip = ip + nExpi
            ip = ip + nExpi
            nHer = ((lb+1)+iAng+2)/2
            nOrder = Max(nHer,nOrder)
            ip = ip + nExpi*3*nHer*(lb+2)
            ip = ip + nExpi*3*nHer*(iAng+1)
            ip = ip + nExpi*3*nHer*(lr+1)
            ip = ip + nExpi*3*nHer*(lb+2)*(iAng+1)*(lr+1)
            ip = ip + nExpi
*
            MmSROG = Max(MmSROG,ip)
            ip = ip - nExpi
     &         * (6 + 3*nHer*((lb+2) + (iAng+1) + (lr+1)
     &         +  (lb+2)*(iAng+1)*(lr+1)) + 1)
*
            ip = ip + Max(nExpi*nac,ncb*nExpi)
            MmSROG = Max(MmSROG,ip)
*
 1966    Continue
 1960 Continue
*
      Return
      End
