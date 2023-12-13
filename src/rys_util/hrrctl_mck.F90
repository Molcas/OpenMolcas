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

subroutine HrrCtl_mck(Arr1,nArr1,Arr2,nArr2,la,lb,lc,ld,nabMax,ncdMax,nTR,A,B,C,D,IfHss,IfGrd)
!***********************************************************************
!                                                                      *
! Object: to act as a shell towards the HRR subroutines.               *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, SWEDEN                               *
!***********************************************************************

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nArr1, nArr2, la, lb, lc, ld, nabMax, ncdMax, nTR
real(kind=wp), intent(inout) :: Arr1(nTR,3*nArr1)
real(kind=wp), intent(out) :: Arr2(nTR,3*nArr2)
real(kind=wp), intent(in) :: A(3), B(3), C(3), D(3)
logical(kind=iwp), intent(in) :: IfHss(4,3,4,3), IfGrd(3,4)

call Hrr2Da_mck(Arr1,nTR,nabMax,ncdMax,Arr2,A,B,la,lb,lc,ld,IfHss,IfGrd)

call Hrr2Db_mck(Arr2,nTR,ncdMax,Arr1,C,D,la,lb,lc,ld,IfHss,IfGrd)

return

end subroutine HrrCtl_mck
