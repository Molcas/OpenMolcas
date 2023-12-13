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
! Copyright (C) 1992,1993, Roland Lindh                                *
!***********************************************************************
!define _DEBUGPRINT_
      SubRoutine Screen(iOffZ,iOffE,nZeta,nEta,mZeta,mEta,lZeta,lEta,
     &                  k2Data1,k2Data2,
     &                  Zeta,ZInv,P,KappAB,IndZet,
     &                  Eta,EInv,Q,KappCD,IndEta,
     &                  Dij,Dkl,
     &                  iphX1,iphY1,iphZ1,iphX2,iphY2,iphZ2,CutDInt,
     &                  CutInt,vij,vkl,vik,vil,vjk,vjl,
     &                  Prescreen_On_Int_Only)
!***********************************************************************
!                                                                      *
! Object: to prescreen the integrals for Direct SCF                    *
!                                                                      *
!   nZeta, nEta : unpartitioned length of primitives.                  *
!                                                                      *
!   mZeta, mEta : section length due to partitioning. These are usually*
!                 equal to nZeta and nEta.                             *
!                                                                      *
!   lZeta, lEta : section length after prescreening.                   *
!                                                                      *
!   Observe that the integrals are ordered for optimal performance in  *
!   the integral generation step whereas the 1st order density is      *
!   canonically ordered.                                               *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, SWEDEN                               *
!             March '92                                                *
!             April '92 modified for gradient estimate                 *
!             January '93 modified to Direct SCF                       *
!***********************************************************************
      use Constants
      use k2_structure, only: k2_type
      Implicit None
      Integer :: iOffZ, iOffE, mZeta, nZeta, mEta, nEta
      Real*8, Intent(out) :: Zeta(mZeta), ZInv(mZeta), KappAB(mZeta),
     &                       P(nZeta,3),
     &                       Eta(mEta),EInv(mEta), KappCD(mEta),
     &                       Q(nEta,3)
      Type(k2_type), intent(in):: k2Data1, k2Data2
      Real*8 Dij(nZeta), Dkl(nEta)
      Integer, Intent(out) :: lZeta, lEta, IndZet(nZeta), IndEta(nEta)
      Integer  iphX1,iphY1,iphZ1,iphX2,iphY2,iphZ2
      Real*8 CutDInt,CutInt,vij,vkl,vik,vil,vjk,vjl
      Logical Prescreen_On_Int_Only
      ![all the others are intent(in)]
!
!     decalaration of local variables...
      Integer iZeta, jZeta, iEta, jEta
      Real*8 abMax,cdMax, DMax, Cut, vMax, ppaa, aaaa, test

      abMax = k2Data1%abMax
      cdMax = k2Data2%abMax
!
      If (.Not.Prescreen_On_Int_Only) Then
!
!        The former technique will be used for integrals which will
!        be partially stored.
!
!-----   Set up prescreening parameters
!
         DMax = Max(vik,vil,vjk,vjl)/Four            ! AO basis
!
!-----   Combine threshold with maximum of density contribution.
!        Note that we are using the square to avoid doing the
!        square root in the denominator!
         Cut  = Max(CutDInt/Max(DMax,vij,vkl),CutDInt)      ! AO basis
!
!-----   Test on the largest integral which will be produced by this
!        batch.
!
         vMax = abMax*cdMax                     ! integral in AO basis
!-----   Skip prescreening if too low.
         If (vMax.lt.Cut) Then
            lZeta = 0
            lEta  = 0
            Go To 999
         End If
      Else
         DMax=Zero
      End If
!
      lZeta=0
      lEta =0
!
!---- prescreen zeta pairs
      If (.Not.Prescreen_On_Int_Only) Then
         lZeta=0
         Do iZeta = 1, mZeta
            ppaa= k2Data1%ab(iOffZ+iZeta) * cdMax
            aaaa= k2Data1%abCon(iOffZ+iZeta) * cdMax
            jZeta = k2Data1%IndZ(iOffZ+iZeta)
            Test=ppaa*Dij(jZeta)+aaaa*(DMax+vkl)
            If (Test.ge.CutDInt) Then
               lZeta=lZeta+1
               Zeta(lZeta)  = k2Data1%Zeta(iOffZ+iZeta)
               KappAB(lZeta)= k2Data1%Kappa(iOffZ+iZeta)
               P(lZeta,1)   = k2Data1%Pcoor(iOffZ+iZeta,1)
               P(lZeta,2)   = k2Data1%Pcoor(iOffZ+iZeta,2)
               P(lZeta,3)   = k2Data1%Pcoor(iOffZ+iZeta,3)
               IndZet(lZeta)= jZeta
               ZInv(lZeta)  = k2Data1%ZInv(iOffZ+iZeta)
            End If
         End Do
      Else
         lZeta=0
         Do iZeta = 1, mZeta
            aaaa= k2Data1%abCon(iOffZ+iZeta) * cdMax
            If (aaaa.ge.CutInt) Then
               lZeta=lZeta+1
               Zeta(lZeta)  = k2Data1%Zeta(iOffZ+iZeta)
               KappAB(lZeta)= k2Data1%Kappa(iOffZ+iZeta)
               P(lZeta,1)   = k2Data1%Pcoor(iOffZ+iZeta,1)
               P(lZeta,2)   = k2Data1%Pcoor(iOffZ+iZeta,2)
               P(lZeta,3)   = k2Data1%Pcoor(iOffZ+iZeta,3)
               IndZet(lZeta)= k2Data1%IndZ(iOffZ+iZeta)
               ZInv(lZeta)  = k2Data1%ZInv(iOffZ+iZeta)
            End If
         End Do
      End If
      If (lZeta.eq.0) Go To 999
!
      If (iphX1.ne.1) Call DScal_(lZeta,-One,P(1,1),1)
      If (iphY1.ne.1) Call DScal_(lZeta,-One,P(1,2),1)
      If (iphZ1.ne.1) Call DScal_(lZeta,-One,P(1,3),1)
!
!---- prescreen eta pairs
      If (.Not.Prescreen_On_Int_Only) Then
         lEta=0
         Do iEta = 1, mEta
            ppaa= k2Data2%ab(iOffE+iEta) * abMax
            aaaa= k2Data2%abCon(iOffE+iEta) * abMax
            jEta = k2Data2%IndZ(iOffE+iEta)
            Test=ppaa*Dkl(jEta)+aaaa*(DMax+vij)
            If (Test.ge.CutDInt) Then
               lEta=lEta+1
               IndEta(lEta)= jEta
               Eta(lEta)   = k2Data2%Zeta(iOffE+iEta)
               KappCD(lEta)= k2Data2%Kappa(iOffE+iEta)
               Q(lEta,1)   = k2Data2%Pcoor(iOffE+iEta,1)
               Q(lEta,2)   = k2Data2%Pcoor(iOffE+iEta,2)
               Q(lEta,3)   = k2Data2%PCoor(iOffE+iEta,3)
               EInv(lEta)  = k2Data2%ZInv(iOffE+iEta)
            End If
         End Do
      Else
         lEta=0
         Do iEta = 1, mEta
            aaaa= k2Data2%abCon(iOffE+iEta) * abMax
            If (aaaa.ge.CutInt) Then
               lEta=lEta+1
               IndEta(lEta)= k2Data2%IndZ(iOffE+iEta)
               Eta(lEta)   = k2Data2%Zeta(iOffE+iEta)
               KappCD(lEta)= k2Data2%Kappa(iOffE+iEta)
               Q(lEta,1)   = k2Data2%Pcoor(iOffE+iEta,1)
               Q(lEta,2)   = k2Data2%Pcoor(iOffE+iEta,2)
               Q(lEta,3)   = k2Data2%PCoor(iOffE+iEta,3)
               EInv(lEta)  = k2Data2%ZInv(iOffE+iEta)
            End If
         End Do
      End If
      If (lEta.eq.0) Go To 999
!
      If (iphX2.ne.1) Call DScal_(lEta,-One,Q(1,1),1)
      If (iphY2.ne.1) Call DScal_(lEta,-One,Q(1,2),1)
      If (iphZ2.ne.1) Call DScal_(lEta,-One,Q(1,3),1)
!
 999  Continue
!define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
      Write (6,*) ' In Screen'
      Call RecPrt(' Zeta  ',' ',Zeta,lZeta,1)
      Call RecPrt(' Eta   ',' ', Eta, lEta,1)
      Call RecPrt(' ZInv  ',' ',ZInv,lZeta,1)
      Call RecPrt(' EInv  ',' ',EInv,lEta, 1)
      Call RecPrt(' Px    ',' ',P(1,1),lZeta,1)
      Call RecPrt(' Py    ',' ',P(1,2),lZeta,1)
      Call RecPrt(' Pz    ',' ',P(1,3),lZeta,1)
      Call RecPrt(' Qx    ',' ',Q(1,1),lEta, 1)
      Call RecPrt(' Qy    ',' ',Q(1,2),lEta, 1)
      Call RecPrt(' Qz    ',' ',Q(1,3),lEta, 1)
      Call RecPrt(' KappAB',' ',KappAB,lZeta,1)
      Call RecPrt(' KappCD',' ',KappCD,lEta, 1)
#endif
      End SubRoutine Screen
