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

use Definitions, only: iwp, u6

implicit none
#include "main_cvb.fh"
integer(kind=iwp) :: nconf1, iconfs(noe,nconf1), nel1
#include "WrkSpc.fh"
integer(kind=iwp) :: i1, iconf, ii, ioffs, iorb
integer(kind=iwp), external :: mstacki_cvb

i1 = mstacki_cvb(noe)
! Main loop over configurations:
do iconf=1,nconf1
  ! Prepare iwork(i1) for print
  ioffs = i1-1
  do iorb=1,norb
    if (iconfs(iorb,iconf) == 2) then
      iwork(1+ioffs) = iorb
      iwork(2+ioffs) = iorb
      ioffs = ioffs+2
    end if
  end do
  do iorb=1,norb
    if (iconfs(iorb,iconf) == 1) then
      iwork(1+ioffs) = iorb
      ioffs = ioffs+1
    end if
  end do
  write(u6,'(i8,a,20i3)') iconf,'   =>  ',(iwork(ii+i1-1),ii=1,nel1)
end do
call mfreei_cvb(i1)

return

end subroutine cnfprt_cvb
