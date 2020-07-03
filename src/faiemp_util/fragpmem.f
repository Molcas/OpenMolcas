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
      Subroutine FragPMem(nHer,MemFrag,la,lb,lr)
************************************************************************
*  Object: to compute the number of real*8 the kernel routine will     *
*          need for the computation of a matrix element between two    *
*          cartesian Gaussin functions with the total angular momentum *
*          of la and lb (la=0 s-function, la=1 p-function, etc.)       *
*          lr is the order of the operator (this is only used when the *
*          integrals are computed with the Hermite-Gauss quadrature).  *
*                                                                      *
*  Called from: OneEl                                                  *
*                                                                      *
*                R. Lindh, Dept. of Theor. Chem. Univ. of Lund, Sweden *
************************************************************************
      Use Basis_Info
      Implicit None
#include "itmax.fh"
#include "info.fh"
      Integer nHer, MemFrag, la, lb, lr,
     &        i, nElem, maxDensSize, iCnttp, jCnttp, iAng, jAng,
     &        iShll, jShll, ip, nac, ncb, MemMlt, nH
* statement function
      nElem(i) = (i+1)*(i+2)/2
*
      MemFrag = 0
      maxDensSize = 0
c        largest possible fragment energy weighted density matrix
      Do iCnttp = 1, nCnttp
        If(dbsc(iCnttp)%nFragType.gt.0) maxDensSize = Max(maxDensSize,
     &                        dbsc(iCnttp)%nFragDens
     &                      *(dbsc(iCnttp)%nFragDens+1)/2)
      End Do
c
      Do iCnttp = 1, nCnttp
      If (.Not.FragCnttp(iCnttp)) cycle
c
         Do iAng = 0, nVal_Shells(iCnttp)-1
         iShll = ipVal(iCnttp) + iAng
         If (nExp(iShll).eq.0 .or. nBasis(iShll).eq.0) cycle
*
            Do jCnttp = iCnttp, nCnttp
* still figure out how to loop only over the centers belonging to the
* same fragment (keep track of mdc?) ! still to be done !!!
            If (.Not.FragCnttp(jCnttp)) cycle
c
               Do jAng = 0, nVal_Shells(jCnttp)-1
               jShll = ipVal(jCnttp) + jAng
               If (nExp(jShll).eq.0 .or. nBasis(jShll).eq.0) cycle
!              Go To 1976
*
               ip = 2 * maxDensSize
              nac = nElem(la)*nElem(iAng)
               ip = ip + nExp(iShll)*nac
               ip = ip + 3 * nExp(iShll)
               ip = ip + nExp(iShll)
               ip = ip + nExp(iShll)
               ip = ip + nExp(iShll)
*
              Call MltMmP(nH,MemMlt,la,iAng,lr)
             nHer = Max(nH,nHer)
          MemFrag = Max(MemFrag,ip+nExp(iShll)*MemMlt)
               ip = ip - 6 * nExp(iShll)
*
              ncb = nElem(jAng)*nElem(lb)
               ip = ip + nExp(jShll)*ncb
               ip = ip + 3 * nExp(jShll)
               ip = ip + nExp(jShll)
               ip = ip + nExp(jShll)
               ip = ip + nExp(jShll)
*
              Call MltMmP(nH,MemMlt,jAng,lb,lr)
             nHer = Max(nH,nHer)
          MemFrag = Max(MemFrag,ip+nExp(jShll)*MemMlt)
               ip = ip - 6 * nExp(jShll)
*
               ip = ip + Max( nac * Max(nExp(iShll),nBasis(jShll)),
     &                        ncb * nBasis(jShll) )
          MemFrag = Max(MemFrag,ip)
*
c 1976          Continue
               Enddo
c
c 1970       Continue
            Enddo
c
c 1966    Continue
         Enddo
c
c 1960 Continue
      Enddo
*
      Return
      End
