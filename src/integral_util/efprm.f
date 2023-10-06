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
! Copyright (C) 1991,1995, Roland Lindh                                *
!***********************************************************************
!#define _DEBUGPRINT_
      SubRoutine EFPrm(Alpha,nAlpha,Beta, nBeta,Zeta,ZInv,rKappa,P,
     &                 Final,nZeta,nComp,la,lb,A,RB,nRys,
     &                 Array,nArr,Ccoor,nOrdOp)
!***********************************************************************
!                                                                      *
! Object: kernel routine for the computation of electric field         *
!         integrals.                                                   *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry, University *
!             of Lund, Sweden, January '91                             *
!                                                                      *
! Modified for explicit code, R. Lindh, February '95.                  *
!***********************************************************************
      use Constants, only: Zero, One
      Implicit None
      External TNAI, Fake, XCff2D, XRys2D
      Integer nZeta, la, lb, nComp, nAlpha, nBeta, nRys, nArr, nOrdOp
      Real*8 Final(nZeta,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2,nComp),
     &       Zeta(nZeta), ZInv(nZeta), Alpha(nAlpha), Beta(nBeta),
     &       rKappa(nZeta), P(nZeta,3), A(3), RB(3),
     &       Array(nZeta*nArr), Ccoor(3)

!---- Local Arrays
      Real*8 Coori(3,4), CoorAC(3,2)
      Logical EQ, NoSpecial
      Integer iAnga(4)
      Integer ixyz, nElem, nabSz
      Integer mabMin, mabMax, mcdMin, mcdMax, lab, kab, lcd, labcd,
     &        ip1, ip2, nMem, mArr, nT, ip3, ipIn, nFlop
#ifdef _DEBUGPRINT_
      Integer iElem, jElem
      Character*80 Label
#endif
!
!     Statement function for Cartesian index
!
      nElem(ixyz) = (ixyz+1)*(ixyz+2)/2
      nabSz(ixyz) = (ixyz+1)*(ixyz+2)*(ixyz+3)/6  - 1
!
#ifdef _DEBUGPRINT_
      Call RecPrt(' In EFPrm: Alpha',' ',Alpha,nAlpha,1)
      Call RecPrt(' In EFPrm: Beta',' ',Beta,nBeta,1)
#else
      If (.False.) Call Unused_Real_Array(Alpha)
      If (.False.) Call Unused_Real_Array(Beta)
#endif
!
      Final(:,:,:,:)=Zero
!
      iAnga(1) = la
      iAnga(2) = lb
      iAnga(3) = nOrdOp
      iAnga(4) = 0
      call dcopy_(3, A,1,Coori(1,1),1)
      call dcopy_(3,RB,1,Coori(1,2),1)
      mabMin=nabSz(Max(la,lb)-1)+1
      mabMax=nabSz(la+lb)
      If (EQ(A,RB)) mabMin=nabSz(la+lb-1)+1
      mcdMin=nabSz(nOrdOp-1)+1
      mcdMax=nabSz(nOrdOp)
      lab=(mabMax-mabMin+1)
      kab=nElem(la)*nElem(lb)
      lcd=(mcdMax-mcdMin+1)
      labcd=lab*lcd
!
!---- Compute Flop's and size of work array which HRR will Use.
!
      Call mHRR(la,lb,nFLOP,nMem)
!
!---- Distribute the work array
!
      ip2 = 1
      ip1 = ip2 + nZeta*Max(labcd,lcd*nMem)
      mArr = nArr - Max(labcd,lcd*nMem)
!
!---- Find center to accumulate angular momentum on. (HRR)
!
      If (la.ge.lb) Then
         call dcopy_(3, A,1,CoorAC(1,1),1)
      Else
         call dcopy_(3,RB,1,CoorAC(1,1),1)
      End If
!
      call dcopy_(3,Ccoor,1,CoorAC(1,2),1)
      call dcopy_(3,Ccoor,1, Coori(1,3),1)
      call dcopy_(3,Ccoor,1, Coori(1,4),1)
!
!-----Compute integrals with the Rys-Gauss quadrature.
!
      nT=nZeta
      NoSpecial=.True.
      Call Rys(iAnga,nT,Zeta,ZInv,nZeta,
     &         [One],[One],1,P,nZeta,
     &         CCoor,1,rKappa,[One],Coori,Coori,CoorAC,
     &         mabMin,mabMax,mcdMin,mcdMax,Array(ip1),mArr*nZeta,
     &         TNAI,Fake,XCff2D,XRys2D,NoSpecial)
!
!---- The integrals are now ordered as ijkl,e,f
!
!     a) Change the order to f,ijkl,e
!     b) Unfold e to ab, f,ijkl,ab
!     c) Change the order back to ijkl,ab,f
!
!a)--
!
      Call DGetMO(Array(ip1),nZeta*lab,nZeta*lab,lcd,Array(ip2),lcd)
!
!b)-- Use the HRR to unfold e to ab
!
      Call HRR(la,lb,A,RB,Array(ip2),lcd*nZeta,nMem,ipIn)
      ip3=ip2-1+ipIn
!
!c)--
!
      Call DGetMO(Array(ip3),lcd,lcd,nZeta*kab,Final,nZeta*kab)
      Call DScal_(nZeta*kab*lcd,-One,Final,1)
!
#ifdef _DEBUGPRINT_
      Write (6,*) ' In EFPrm la,lb=',la,lb
      Do 400 iElem = 1, nElem(la)
         Do 410 jElem = 1, nElem(lb)
            If (lcd.eq.1) Then
               Write (Label,'(A,I2,A,I2,A)')
     &               ' EFPrm: Final (',iElem,',',jElem,') '
               Call RecPrt(Label,' ',Final(1,iElem,jElem,1),nZeta,1)
            Else If (lcd.eq.3) tHEN
               Write (Label,'(A,I2,A,I2,A)')
     &               ' EFPrm: Final (',iElem,',',jElem,',x) '
               Call RecPrt(Label,' ',Final(1,iElem,jElem,1),nZeta,1)
               Write (Label,'(A,I2,A,I2,A)')
     &               ' EFPrm: Final (',iElem,',',jElem,',y) '
               Call RecPrt(Label,' ',Final(1,iElem,jElem,2),nZeta,1)
               Write (Label,'(A,I2,A,I2,A)')
     &               ' EFPrm: Final (',iElem,',',jElem,',z) '
               Call RecPrt(Label,' ',Final(1,iElem,jElem,3),nZeta,1)
            End If
 410     Continue
 400  Continue
#endif
!
      Return
! Avoid unused argument warnings
      If (.False.) Call Unused_integer(nRys)
      End SubRoutine EFPrm
