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

subroutine Inter(iSet1,nSet1,iSet2,nSet2,iInter,nInter)
!***********************************************************************
!                                                                      *
! Object : to form the intersection of two sets.                       *
!                                                                      *
! Called from: NucAtt                                                  *
!              TwoEl                                                   *
!                                                                      *
! Calling    : None                                                    *
!                                                                      *
!     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
!             February '90                                             *
!***********************************************************************

use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: nSet1, iSet1(0:nSet1-1), nSet2, iSet2(0:nSet2-1)
integer(kind=iwp), intent(inout) :: iInter(0:7)
integer(kind=iwp), intent(out) :: nInter
integer(kind=iwp) :: i1, i2

nInter = 0
do i1=0,nSet1-1
  do i2=0,nSet2-1
    if (iSet1(i1) == iSet2(i2)) then
      iInter(nInter) = iSet1(i1)
      nInter = nInter+1
      exit
    end if
  end do
end do

return

end subroutine Inter
