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
! Copyright (C) 2006, Thomas Bondo Pedersen                            *
!***********************************************************************

subroutine Domain_Histo1(iDomain,nAtom,nOcc,iCount,i_min,i_max)

use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: nAtom, nOcc, iDomain(0:nAtom,nOcc), i_min, i_max
integer(kind=iwp), intent(out) :: iCount(i_max-i_min+1)
integer(kind=iwp) :: i, iC, nC

nC = i_max-i_min+1
iCount(1:nC) = 0

do i=1,nOcc
  iC = iDomain(0,i)-i_min+1
  iCount(iC) = iCount(iC)+1
end do

end subroutine Domain_Histo1
