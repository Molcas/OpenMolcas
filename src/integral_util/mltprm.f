!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 1990, Roland Lindh                                     *
!               1990, IBM                                              *
!***********************************************************************
!#define _DEBUGPRINT_
      SubRoutine MltPrm(Alpha,nAlpha,Beta, nBeta,Zeta,ZInv,rKappa,P,
     &                  Final,nZeta,nComp,la,lb,A,RB,nHer,
     &                  Array,nArr,Ccoor,nOrdOp)
!***********************************************************************
!                                                                      *
! Object: to compute the multipole moments integrals with the          *
!         Gauss-Hermite quadrature.                                    *
!                                                                      *
!     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
!             November '90                                             *
!             Modified to multipole moments November '90               *
!***********************************************************************
      use Her_RW, only: HerR, HerW, iHerR, iHerw
      use Constants
      Implicit None
      Integer nZeta, la, lb, nComp, nAlpha, nBeta, nArr, nHer, nOrdOp
      Real*8 Final(nZeta,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2,nComp),
     &       Zeta(nZeta), ZInv(nZeta), Alpha(nAlpha), Beta(nBeta),
     &       rKappa(nZeta), P(nZeta,3), A(3), RB(3),
     &       Array(nZeta*nArr), Ccoor(3)

      Logical ABeq(3)
      Integer nip, ipAxyz, ipBxyz, ipRxyz, ipQxyz
!
!     Statement function for Cartesian index
!
!     nElem(i) = (i+1)*(i+2)/2
!     Ind(ixyz,ix,iz) = (ixyz-ix)*(ixyz-ix+1)/2 + iz + 1
!     iOff(ixyz) = ixyz*(ixyz+1)*(ixyz+2)/6
!     Index(ixyz,ix,iz) = Ind(ixyz,ix,iz) + iOff(ixyz)
!
      ABeq(:) = A(:).eq.RB(:)
!
      nip = 1
      ipAxyz = nip
      nip = nip + nZeta*3*nHer*(la+1)
      ipBxyz = nip
      nip = nip + nZeta*3*nHer*(lb+1)
      ipRxyz = nip
      nip = nip + nZeta*3*nHer*(nOrdOp+1)
      ipQxyz = nip
      nip = nip + nZeta*3*(la+1)*(lb+1)*(nOrdOp+1)
      If (nip-1.gt.nArr*nZeta) Then
         Call WarningMessage(2,'MltPrm: nip-1.gt.nArr*nZeta')
         Write (6,*) ' nArr is Wrong! ', nip-1,' > ',nArr*nZeta
         Write (6,*) ' Abend in MltPrm'
         Call Abend()
      End If
!
#ifdef _DEBUGPRINT_
      Call RecPrt(' In MltPrm: A',' ',A,1,3)
      Call RecPrt(' In MltPrm: RB',' ',RB,1,3)
      Call RecPrt(' In MltPrm: Ccoor',' ',Ccoor,1,3)
      Call RecPrt(' In MltPrm: Kappa',' ',rKappa,nAlpha,nBeta)
      Call RecPrt(' In MltPrm: Zeta',' ',Zeta,nAlpha,nBeta)
      Call RecPrt(' In MltPrm: P',' ',P,nZeta,3)
      Write (6,*) ' In MltPrm: la,lb=',la,lb
#endif
!
!     Compute the cartesian values of the basis functions angular part
!
      Call CrtCmp(Zeta,P,nZeta,A,Array(ipAxyz),
     &               la,HerR(iHerR(nHer)),nHer,ABeq)
      Call CrtCmp(Zeta,P,nZeta,RB,Array(ipBxyz),
     &               lb,HerR(iHerR(nHer)),nHer,ABeq)
!
!     Compute the contribution from the multipole moment operator
!
      ABeq(:) = .False.
      Call CrtCmp(Zeta,P,nZeta,Ccoor,Array(ipRxyz),
     &            nOrdOp,HerR(iHerR(nHer)),nHer,ABeq)
!
!     Compute the cartesian components for the multipole moment
!     integrals. The integrals are factorized into components.
!
       Call Assmbl(Array(ipQxyz),
     &             Array(ipAxyz),la,
     &             Array(ipRxyz),nOrdOp,
     &             Array(ipBxyz),lb,
     &             nZeta,HerW(iHerW(nHer)),nHer)
!
!     Combine the cartesian components to the full one electron
!     integral.
!
      Call CmbnMP(Array(ipQxyz),nZeta,la,lb,nOrdOp,Zeta,rKappa,
     &            Final,nComp)
!
      Return
! Avoid unused argument warnings
      If (.False.) Then
         Call Unused_real_array(Alpha)
         Call Unused_real_array(Beta)
         Call Unused_real_array(ZInv)
      End If
      End SubRoutine MltPrm
