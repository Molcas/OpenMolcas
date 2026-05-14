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

!#define _DEBUGPRINT_
subroutine Screen(iOffZ,iOffE,nZeta,nEta,mZeta,mEta,lZeta,lEta,k2Data1,k2Data2,Zeta,ZInv,P,KappAB,IndZet,Eta,EInv,Q,KappCD,IndEta, &
                  Dij,Dkl,iphX1,iphY1,iphZ1,iphX2,iphY2,iphZ2,CutDInt,CutInt,vij,vkl,vik,vil,vjk,vjl,Prescreen_On_Int_Only)
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

use k2_structure, only: k2_type
use Constants, only: Zero, Four
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: iOffZ, iOffE, nZeta, nEta, mZeta, mEta, iphX1, iphY1, iphZ1, iphX2, iphY2, iphZ2
integer(kind=iwp), intent(out) :: lZeta, lEta, IndZet(nZeta), IndEta(nEta)
type(k2_type), intent(in) :: k2Data1, k2Data2
real(kind=wp), intent(out) :: Zeta(mZeta), ZInv(mZeta), P(nZeta,3), KappAB(mZeta), Eta(mEta), EInv(mEta), Q(nEta,3), KappCD(mEta)
real(kind=wp), intent(in) :: Dij(nZeta), Dkl(nEta), CutDInt, CutInt, vij, vkl, vik, vil, vjk, vjl
logical(kind=iwp), intent(in) :: Prescreen_On_Int_Only
integer(kind=iwp) :: iEta, iZeta, jEta, jZeta
real(kind=wp) :: aaaa, abMax, cdMax, Cut, DMax, ppaa, test, vMax
logical(kind=iwp) :: Skip

abMax = k2Data1%abConMax
cdMax = k2Data2%abConMax

Skip = .false.
if (.not. Prescreen_On_Int_Only) then

  ! The former technique will be used for integrals which will
  ! be partially stored.

  ! Set up prescreening parameters

  DMax = max(vik,vil,vjk,vjl)/Four            ! AO basis

  ! Combine threshold with maximum of density contribution.
  ! Note that we are using the square to avoid doing the
  ! square root in the denominator!
  Cut = max(CutDInt/max(DMax,vij,vkl),CutDInt)      ! AO basis

  ! Test on the largest integral which will be produced by this batch.

  vMax = abMax*cdMax                     ! integral in AO basis
  ! Skip prescreening if too low.
  if (vMax < Cut) Skip = .true.
else
  DMax = Zero
end if

lZeta = 0
lEta = 0

if (.not. Skip) then
  ! prescreen zeta pairs
  if (.not. Prescreen_On_Int_Only) then
    do iZeta=1,mZeta
      ppaa = k2Data1%ab(iOffZ+iZeta)*cdMax
      aaaa = k2Data1%abCon(iOffZ+iZeta)*cdMax
      jZeta = k2Data1%IndZ(iOffZ+iZeta)
      Test = ppaa*Dij(jZeta)+aaaa*(DMax+vkl)
      if (Test >= CutDInt) then
        lZeta = lZeta+1
        Zeta(lZeta) = k2Data1%Zeta(iOffZ+iZeta)
        KappAB(lZeta) = k2Data1%Kappa(iOffZ+iZeta)
        P(lZeta,1) = k2Data1%Pcoor(iOffZ+iZeta,1)
        P(lZeta,2) = k2Data1%Pcoor(iOffZ+iZeta,2)
        P(lZeta,3) = k2Data1%Pcoor(iOffZ+iZeta,3)
        IndZet(lZeta) = jZeta
        ZInv(lZeta) = k2Data1%ZInv(iOffZ+iZeta)
      end if
    end do
  else
    do iZeta=1,mZeta
      aaaa = k2Data1%abCon(iOffZ+iZeta)*cdMax
      if (aaaa >= CutInt) then
        lZeta = lZeta+1
        Zeta(lZeta) = k2Data1%Zeta(iOffZ+iZeta)
        KappAB(lZeta) = k2Data1%Kappa(iOffZ+iZeta)
        P(lZeta,1) = k2Data1%Pcoor(iOffZ+iZeta,1)
        P(lZeta,2) = k2Data1%Pcoor(iOffZ+iZeta,2)
        P(lZeta,3) = k2Data1%Pcoor(iOffZ+iZeta,3)
        IndZet(lZeta) = k2Data1%IndZ(iOffZ+iZeta)
        ZInv(lZeta) = k2Data1%ZInv(iOffZ+iZeta)
      end if
    end do
  end if

  if (lZeta /= 0) then
    if (iphX1 /= 1) P(1:lZeta,1) = -P(1:lZeta,1)
    if (iphY1 /= 1) P(1:lZeta,2) = -P(1:lZeta,2)
    if (iphZ1 /= 1) P(1:lZeta,3) = -P(1:lZeta,3)

    ! prescreen eta pairs
    if (.not. Prescreen_On_Int_Only) then
      do iEta=1,mEta
        ppaa = k2Data2%ab(iOffE+iEta)*abMax
        aaaa = k2Data2%abCon(iOffE+iEta)*abMax
        jEta = k2Data2%IndZ(iOffE+iEta)
        Test = ppaa*Dkl(jEta)+aaaa*(DMax+vij)
        if (Test >= CutDInt) then
          lEta = lEta+1
          IndEta(lEta) = jEta
          Eta(lEta) = k2Data2%Zeta(iOffE+iEta)
          KappCD(lEta) = k2Data2%Kappa(iOffE+iEta)
          Q(lEta,1) = k2Data2%Pcoor(iOffE+iEta,1)
          Q(lEta,2) = k2Data2%Pcoor(iOffE+iEta,2)
          Q(lEta,3) = k2Data2%PCoor(iOffE+iEta,3)
          EInv(lEta) = k2Data2%ZInv(iOffE+iEta)
        end if
      end do
    else
      do iEta=1,mEta
        aaaa = k2Data2%abCon(iOffE+iEta)*abMax
        if (aaaa >= CutInt) then
          lEta = lEta+1
          IndEta(lEta) = k2Data2%IndZ(iOffE+iEta)
          Eta(lEta) = k2Data2%Zeta(iOffE+iEta)
          KappCD(lEta) = k2Data2%Kappa(iOffE+iEta)
          Q(lEta,1) = k2Data2%Pcoor(iOffE+iEta,1)
          Q(lEta,2) = k2Data2%Pcoor(iOffE+iEta,2)
          Q(lEta,3) = k2Data2%PCoor(iOffE+iEta,3)
          EInv(lEta) = k2Data2%ZInv(iOffE+iEta)
        end if
      end do
    end if
    if (lEta /= 0) then
      if (iphX2 /= 1) Q(1:lEta,1) = -Q(1:lEta,1)
      if (iphY2 /= 1) Q(1:lEta,2) = -Q(1:lEta,2)
      if (iphZ2 /= 1) Q(1:lEta,3) = -Q(1:lEta,3)
    end if
  end if
end if

!define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
write(u6,*) ' In Screen'
call RecPrt(' Zeta  ',' ',Zeta,lZeta,1)
call RecPrt(' Eta   ',' ',Eta,lEta,1)
call RecPrt(' ZInv  ',' ',ZInv,lZeta,1)
call RecPrt(' EInv  ',' ',EInv,lEta,1)
call RecPrt(' Px    ',' ',P(1,1),lZeta,1)
call RecPrt(' Py    ',' ',P(1,2),lZeta,1)
call RecPrt(' Pz    ',' ',P(1,3),lZeta,1)
call RecPrt(' Qx    ',' ',Q(1,1),lEta,1)
call RecPrt(' Qy    ',' ',Q(1,2),lEta,1)
call RecPrt(' Qz    ',' ',Q(1,3),lEta,1)
call RecPrt(' KappAB',' ',KappAB,lZeta,1)
call RecPrt(' KappCD',' ',KappCD,lEta,1)
#endif

end subroutine Screen
