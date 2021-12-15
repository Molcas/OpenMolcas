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
      Integer nHer, MemFrag, la, lb, lr, nExpi, nExpj,
     &        i, nElem, maxDensSize, iCnttp, jCnttp, iAng, jAng,
     &        iShll, jShll, ip, nac, ncb, MemMlt, nH,
     &        nBasisi, nBasisj
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
      If (.Not.dbsc(iCnttp)%Frag) cycle
c
         Do iAng = 0, dbsc(iCnttp)%nVal-1
         iShll = dbsc(iCnttp)%iVal + iAng
         nExpi=Shells(iShll)%nExp
         nBasisi=Shells(iShll)%nBasis
         If (nExpi.eq.0 .or. nBasisi.eq.0) cycle
*
            Do jCnttp = iCnttp, nCnttp
* still figure out how to loop only over the centers belonging to the
* same fragment (keep track of mdc?) ! still to be done !!!
            If (.Not.dbsc(jCnttp)%Frag) cycle
c
               Do jAng = 0, dbsc(jCnttp)%nVal-1
               jShll = dbsc(jCnttp)%iVal + jAng
               nExpj=Shells(jShll)%nExp
               nBasisj=Shells(jShll)%nBasis
               If (nExpj.eq.0 .or. nBasisj.eq.0) cycle
!              Go To 1976
*
               ip = 2 * maxDensSize
              nac = nElem(la)*nElem(iAng)
               ip = ip + nExpi*nac
               ip = ip + 3 * nExpi
               ip = ip + nExpi
               ip = ip + nExpi
               ip = ip + nExpi
*
              Call MltMmP(nH,MemMlt,la,iAng,lr)
             nHer = Max(nH,nHer)
          MemFrag = Max(MemFrag,ip+nExpi*MemMlt)
               ip = ip - 6 * nExpi
*
              ncb = nElem(jAng)*nElem(lb)
               ip = ip + nExpj*ncb
               ip = ip + 3 * nExpj
               ip = ip + nExpj
               ip = ip + nExpj
               ip = ip + nExpj
*
              Call MltMmP(nH,MemMlt,jAng,lb,lr)
             nHer = Max(nH,nHer)
          MemFrag = Max(MemFrag,ip+nExpj*MemMlt)
               ip = ip - 6 * nExpj
*
               ip = ip + Max( nac * Max(nExpi,nBasisj),
     &                        ncb * nBasisj )
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
