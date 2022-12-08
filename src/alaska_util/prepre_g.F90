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

subroutine PrePre_g(nZeta,nEta,mZeta,mEta,lZeta,lEta,Data1,Data2,PreScr,CutGrd)
!***********************************************************************
!                                                                      *
! Object: to preprescreen the integral derivatives.                    *
!                                                                      *
!   nZeta, nEta : unpartitioned length of primitives.                  *
!                                                                      *
!   mZeta, mEta : section length due to partioning. These are usually  *
!                 equal to nZeta and nEta.                             *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, SWEDEN                               *
!             July '92.                                                *
!***********************************************************************

use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nZeta, nEta, mZeta, mEta
integer(kind=iwp), intent(out) :: lZeta, lEta
real(kind=wp), intent(in) :: Data1(nZeta,8), Data2(nEta,8), CutGrd
logical(kind=iwp), intent(out) :: PreScr
integer(kind=iwp) :: iEta, iPrint, iRout, iZeta
real(kind=wp) :: EtaMn, EtaMx, PreMax, PreMin, rKabMn, rKabMx, rKcdMn, rKcdMx, ZetaMn, ZetaMx
#include "print.fh"

iRout = 180
iPrint = nPrint(iRout)
!iQ = 0
if (iPrint >= 99) then
  call RecPrt(' Data1',' ',Data1,nZeta,8)
  call RecPrt(' Data2',' ',Data2,nEta,8)
end if

! Preprescanning

lZeta = mZeta
lEta = mEta
rKabMx = Zero
ZetaMx = Zero
rKabMn = 1.0e72_wp
ZetaMn = Zero
do iZeta=1,mZeta
  if (Data1(iZeta,2) > rKabMx) then
    rKabMx = Data1(iZeta,2)
    ZetaMx = Data1(iZeta,1)
  end if
  if (Data1(iZeta,2) < rKabMn) then
    rKabMn = Data1(iZeta,2)
    ZetaMn = Data1(iZeta,1)
  end if
end do
rKcdMx = Zero
EtaMx = Zero
rKcdMn = 1.0e72_wp
EtaMn = Zero
do iEta=1,mEta
  if (Data2(iEta,2) > rKcdMx) then
    rKcdMx = Data2(iEta,2)
    EtaMx = Data2(iEta,1)
  end if
  if (Data2(iEta,2) < rKcdMn) then
    rKcdMn = Data2(iEta,2)
    EtaMn = Data2(iEta,1)
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
