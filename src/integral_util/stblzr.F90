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

subroutine Stblzr(iU,nU,iV,nV,iR,iM,nM)
!***********************************************************************
!                                                                      *
!  Object: to form the proper stabilizer for a pair of entities        *
!          given the respective stabilizers and the operator acting    *
!          on the second entity.                                       *
!                                                                      *
!     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
!             June '90                                                 *
!***********************************************************************

use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: nU, iU(nU), nV, iV(nV), iR
integer(kind=iwp), intent(out) :: iM(8), nM
integer(kind=iwp) :: i, iRU
logical(kind=iwp) :: UeqV
logical(kind=iwp), external :: RinT

! See if U and V are the same

UeqV = .true.
do i=1,nV
  if (.not. RinT(iU,nU,iV(i))) then
    UeqV = .false.
    exit
  end if
end do
if (UeqV) then
  do i=1,nU
    if (.not. RinT(iV,nV,iU(i))) then
      UeqV = .false.
      exit
    end if
  end do
end if

if (UeqV) then

  ! M is formed as U union RU

  iM(1:nU) = iU(:)
  nM = nU
  do i=1,nU
    iRU = ieor(iR,iU(i))
    if (.not. RinT(iM,nM,iRU)) then
      nM = nM+1
      iM(nM) = iRU
    end if
  end do
else

  ! M is formed as U intersection V

  call Inter(iU,nU,iV,nV,iM,nM)
end if

return

end subroutine Stblzr
