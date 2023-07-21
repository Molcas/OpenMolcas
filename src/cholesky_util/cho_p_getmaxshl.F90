!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

subroutine Cho_P_GetMaxShl(DiaSh,Smax,iShlAB)
!
! Purpose: get max shell pair and update DiaSh.

implicit none
real*8 DiaSh(*)
real*8 Smax
integer iShlAB
#include "cho_para_info.fh"

if (Cho_Real_Par) then
  ! Swap local and global reduced set index arrays and use original
  ! serial routine to get max shell pair.
  call Cho_P_IndxSwp()
  call Cho_GetMaxShl(DiaSh,Smax,iShlAB)
  call Cho_P_IndxSwp()
else
  call Cho_GetMaxShl(DiaSh,Smax,iShlAB)
end if

end subroutine Cho_P_GetMaxShl
