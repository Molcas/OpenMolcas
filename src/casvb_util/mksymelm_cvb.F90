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

subroutine mksymelm_cvb()

use Definitions, only: iwp, u6

implicit none
#include "main_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"
#include "WrkSpc.fh"
integer(kind=iwp) :: ioffs, ishift, isyme
logical(kind=iwp), external :: up2date_cvb ! ... Make: up to date? ...

call rdioff_cvb(8,recinp,ioffs)
call rdr_cvb(work(ls(1)),nsyme*norb*norb,recinp,ioffs)
if ((ip(2) >= 1) .and. (.not. up2date_cvb('PRSYMELM'))) then
  do isyme=1,nsyme
    write(u6,'(/,a,i4,3x,a)') ' Symmetry element no.',isyme,tags(isyme)
    ishift = norb*norb*(isyme-1)
    call mxprint_cvb(work(ishift+ls(1)),norb,norb,0)
  end do
  if (nsyme > 0) write(u6,*) ' '
  call untouch_cvb('PRSYMELM')
end if

return

end subroutine mksymelm_cvb
