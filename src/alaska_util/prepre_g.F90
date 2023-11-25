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
! Copyright (C) 1992, Roland Lindh                                     *
!***********************************************************************

!#define _DEBUGPRINT_
subroutine PrePre_g(mZeta,mEta,lZeta,lEta,PreScr,CutGrd,iOffZ,iOffE,k2Data1,k2Data2)
!***********************************************************************
!                                                                      *
! Object: to preprescreen the integral derivatives.                    *
!                                                                      *
!   mZeta, mEta : section length due to partioning. These are usually  *
!                 equal to nZeta and nEta.                             *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, SWEDEN                               *
!             July '92.                                                *
!***********************************************************************

use k2_structure, only: k2_type
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: mZeta, mEta, iOffZ, iOffE
integer(kind=iwp), intent(out) :: lZeta, lEta
logical(kind=iwp), intent(out) :: PreScr
real(kind=wp), intent(in) :: CutGrd
type(k2_type), intent(in) :: k2Data1, k2Data2
integer(kind=iwp) :: iEta, iZeta
real(kind=wp) :: EtaMn, EtaMx, PreMax, PreMin, rKabMn, rKabMx, rKcdMn, rKcdMx, ZetaMn, ZetaMx

! Preprescanning

lZeta = mZeta
lEta = mEta
rKabMx = Zero
ZetaMx = Zero
rKabMn = 1.0e72_wp
ZetaMn = Zero
do iZeta=1,mZeta
  if (k2Data1%Kappa(iOffZ+iZeta) > rKabMx) then
    rKabMx = k2Data1%Kappa(iOffZ+iZeta)
    ZetaMx = k2Data1%Zeta(iOffZ+iZeta)
  end if
  if (k2Data1%Kappa(iOffZ+iZeta) < rKabMn) then
    rKabMn = k2Data1%Kappa(iOffZ+iZeta)
    ZetaMn = k2Data1%Zeta(iOffZ+iZeta)
  end if
end do
rKcdMx = Zero
EtaMx = Zero
rKcdMn = 1.0e72_wp
EtaMn = Zero
do iEta=1,mEta
  if (k2Data2%Kappa(iOffE+iEta) > rKcdMx) then
    rKcdMx = k2Data2%Kappa(iOffE+iEta)
    EtaMx = k2Data2%Zeta(iOffE+iEta)
  end if
  if (k2Data2%Kappa(iOffE+iEta) < rKcdMn) then
    rKcdMn = k2Data2%Kappa(iOffE+iEta)
    EtaMn = k2Data2%Zeta(iOffE+iEta)
  end if
end do
PreScr = .true.
PreMax = rKabMx*rKcdMx*sqrt(One/(ZetaMx+EtaMx))
PreMin = rKabMn*rKcdMn*sqrt(One/(ZetaMn+EtaMn))
if (PreMin > CutGrd) PreScr = .false.
if (PreMax < 1.0e-4_wp*CutGrd) then
  lZeta = 0
  lEta = 0
end if

return

end subroutine PrePre_g
