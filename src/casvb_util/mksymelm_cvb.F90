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

use casvb_global, only: ipr, norb, nsyme, recinp, symelm, tags
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp) :: ioffs, isyme
logical(kind=iwp), external :: up2date_cvb ! ... Make: up to date? ...

call rdioff_cvb(8,recinp,ioffs)
call rdlow_cvb(symelm,nsyme*norb*norb,recinp,ioffs)
if ((ipr(2) >= 1) .and. (.not. up2date_cvb('PRSYMELM'))) then
  do isyme=1,nsyme
    write(u6,'(/,a,i4,3x,a)') ' Symmetry element no.',isyme,tags(isyme)
    call mxprint_cvb(symelm(:,:,isyme),norb,norb,0)
  end do
  if (nsyme > 0) write(u6,*) ' '
  call untouch_cvb('PRSYMELM')
end if

return

end subroutine mksymelm_cvb
