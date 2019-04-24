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
* Copyright (C) Ben Swerts                                             *
*               2016, Liviu Ungur                                      *
************************************************************************
      Subroutine FragPMmG(nHer,MemFrag,la,lb,lr)
************************************************************************
*  Object: to compute the number of real*8 the kernal routine will     *
*          need for the computation of a matrix element between two    *
*          cartesian Gaussin functions with the total angular momentum *
*          of la and lb (la=0 s-function, la=1 p-function, etc.)       *
*          lr is the order of the operator (this is only used when the *
*          integrals are computed with the Hermite-Gauss quadrature).  *
*                                                                      *
*  Called from: OneEl_g                                                *
*                                                                      *
*  author: Ben Swerts                                                  *
*  modified: Liviu Ungur                                               *
*  based on PrjMmG                                                     *
*                                                                      *
************************************************************************
      Implicit None
#include "itmax.fh"
#include "info.fh"
      Integer nHer,MemFrag,la,lb,lr
      Integer i,nElem,nOrder,maxDensSize
      Integer iCnttp,jCnttp,iAng,jAng,iShll,jShll
      Integer ip,nac,ncb
*
      nElem(i) = (i+1)*(i+2)/2
*
      nOrder = 0
      MemFrag = 0
      maxDensSize = 0
c  largest possible fragment energy weighted density matrix
      Do iCnttp = 1, nCnttp
        If(nFragType(iCnttp).gt.0)
     &  maxDensSize = Max( maxDensSize,
     &                     nFragDens(iCnttp)*(nFragDens(iCnttp)+1)/2 )
      End Do
c
      Do iCnttp = 1, nCnttp
      If (.Not.FragCnttp(iCnttp)) cycle  ! Go To 1960
*
         Do iAng = 0, nVal_Shells(iCnttp)-1
         iShll = ipVal(iCnttp) + iAng
         If (nExp(iShll).eq.0 .or. nBasis(iShll).eq.0) cycle !Go To 1966
*
            Do jCnttp = iCnttp, nCnttp
* still figure out how to loop only over the centers belonging to the
* same fragment (keep track of mdc?) ! still to be done !!!
            If (.Not.FragCnttp(jCnttp)) cycle !Go To 1970
*
               Do jAng = 0, nVal_Shells(jCnttp)-1
               jShll = ipVal(jCnttp) + jAng
               If (nExp(jShll).eq.0 .or. nBasis(jShll).eq.0) cycle
!              Go To 1976
*
               ip =  2 * maxDensSize
              nac =  4 * nElem(la) * nElem(iAng)
               ip = ip + nExp(iShll) * nac
               ip = ip + 3 * nExp(iShll)
               ip = ip + nExp(iShll)
               ip = ip + nExp(iShll)
               ip = ip + nExp(iShll)
             nHer = ((la+1)+iAng+2)/2
           nOrder = Max(nHer,nOrder)
               ip = ip + nExp(iShll) * 3 * nHer * (la+2)
               ip = ip + nExp(iShll) * 3 * nHer * (iAng+1)
               ip = ip + nExp(iShll) * 3 * nHer * (lr+1)
               ip = ip + nExp(iShll) * 3 * nHer * (la+2)*(iAng+1)*(lr+1)
               ip = ip + nExp(iShll)
*
          MemFrag = Max(MemFrag,ip)
               ip = ip - nExp(iShll)
     &            * (6 + 3*nHer*((la+2) + (iAng+1) + (lr+1)
     &            + (la+2)*(iAng+1)*(lr+1)) + 1)
*
              ncb = 4*nElem(jAng)*nElem(lb)
               ip = ip + nExp(jShll)*ncb
               ip = ip + 3 * nExp(jShll)
               ip = ip + nExp(jShll)
               ip = ip + nExp(jShll)
               ip = ip + nExp(jShll)
             nHer = ((lb+1)+jAng+2)/2
           nOrder = Max(nHer,nOrder)
               ip = ip + nExp(jShll)*3*nHer*(lb+2)
               ip = ip + nExp(jShll)*3*nHer*(jAng+1)
               ip = ip + nExp(jShll)*3*nHer*(lr+1)
               ip = ip + nExp(jShll)*3*nHer*(lb+2)*(jAng+1)*(lr+1)
               ip = ip + nExp(jShll)
*
          MemFrag = Max(MemFrag,ip)
               ip = ip - nExp(jShll)
     &            * (6 + 3*nHer*((lb+2) + (jAng+1) + (lr+1)
     &            +  (lb+2)*(jAng+1)*(lr+1)) + 1)
*
               ip = ip + Max(Max(nExp(iShll),nBasis(jShll))*nac,
     &                      ncb*nBasis(jShll))
          MemFrag = Max(MemFrag,ip)
*
c 1976          Continue
               Enddo  !jAng

c 1970       Continue
            Enddo  !jCnttp

c 1966    Continue
         Enddo !iAng

c 1960 Continue
      Enddo ! iCnttp
      nHer = nOrder
*
      Return
      End
