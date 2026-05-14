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
! Copyright (C) 1996-2006, Thorstein Thorsteinsson                     *
!               1996-2006, David L. Cooper                             *
!***********************************************************************

!***********************************************************************
!*                                                                     *
!*  CNFPRT   := Print configurations.                                  *
!*                                                                     *
!***********************************************************************
subroutine cnfprt_cvb(iconfs,nconf1,nel1)

use casvb_global, only: noe, norb
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nconf1, iconfs(noe,nconf1), nel1
integer(kind=iwp) :: iconf, ii, ioffs, iorb
integer(kind=iwp), allocatable :: tmp(:)

call mma_allocate(tmp,noe,label='tmp')
! Main loop over configurations:
do iconf=1,nconf1
  ! Prepare tmp for print
  ioffs = 0
  do iorb=1,norb
    if (iconfs(iorb,iconf) == 2) then
      tmp(1+ioffs) = iorb
      tmp(2+ioffs) = iorb
      ioffs = ioffs+2
    end if
  end do
  do iorb=1,norb
    if (iconfs(iorb,iconf) == 1) then
      tmp(1+ioffs) = iorb
      ioffs = ioffs+1
    end if
  end do
  write(u6,'(i8,a,20i3)') iconf,'   =>  ',(tmp(ii),ii=1,nel1)
end do
call mma_deallocate(tmp)

return

end subroutine cnfprt_cvb
