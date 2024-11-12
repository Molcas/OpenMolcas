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

subroutine Phase(iCmp,jCmp,kCmp,lCmp,iAng,iShll,kOp,ijkl,AOInt)
!***********************************************************************
!                                                                      *
!  Object: To change the phase of the integrals in accordance with the *
!          swapping of the operators operating on the integrals.       *
!                                                                      *
!     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
!             June '90                                                 *
!***********************************************************************

use Index_Functions, only: nTri3_Elem
use Basis_Info, only: Shells
use Real_Spherical, only: iSphCr
use Symmetry_Info, only: iChBas
use Constants, only: One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: iCmp, jCmp, kCmp, lCmp, iAng(4), iShll(4), kOp, ijkl
real(kind=wp), intent(inout) :: AOInt(ijkl,iCmp,jCmp,kCmp,lCmp)
integer(kind=iwp) :: i1, i2, i3, i4, iChBs, ii, jChBs, jj, kChBs, kk, lChBs, ll
real(kind=wp) :: Factor, pa1T, pa2T, pb1T, pb2T
integer(kind=iwp), external :: iPrmt

!call RecPrt(' In Phase: AOInt ',' ',AOInt,ijkl,ijCmp*ijCmp)

! Change phase factor. This is only necessary if T=/=E.

if ((kOp == 0) .or. (iCmp*jCmp*kCmp*lCmp == 0)) return
ii = nTri3_Elem(iAng(1))
jj = nTri3_Elem(iAng(2))
kk = nTri3_Elem(iAng(3))
ll = nTri3_Elem(iAng(4))
do i1=1,iCmp
  iChBs = iChBas(ii+i1)
  if (Shells(iShll(1))%Transf) iChBs = iChBas(iSphCr(ii+i1))
  pa1T = real(iPrmt(kOp,iChBs),kind=wp)
  do i2=1,jCmp
    jChBs = iChBas(jj+i2)
    if (Shells(iShll(2))%Transf) jChBs = iChBas(iSphCr(jj+i2))
    pb1T = real(iPrmt(kOp,jChBs),kind=wp)

    do i3=1,kCmp
      kChBs = iChBas(kk+i3)
      if (Shells(iShll(3))%Transf) kChBs = iChBas(iSphCr(kk+i3))
      pa2T = real(iPrmt(kOp,kChBs),kind=wp)
      do i4=1,lCmp
        lChBs = iChBas(ll+i4)
        if (Shells(iShll(4))%Transf) lChBs = iChBas(iSphCr(ll+i4))
        pb2T = real(iPrmt(kOp,lChBs),kind=wp)
        Factor = pa1T*pb1T*pa2T*pb2T
        if (Factor /= One) AOInt(:,i1,i2,i3,i4) = Factor*AOInt(:,i1,i2,i3,i4)
      end do
    end do
  end do
end do

!call RecPrt(' Exit Phase: AOInt ',' ',AOInt,ijkl,iCmp*jCmp*kCmp*lCmp)

return

end subroutine Phase
