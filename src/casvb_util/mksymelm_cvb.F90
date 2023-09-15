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

implicit real*8(a-h,o-z)
! ... Make: up to date? ...
logical, external :: up2date_cvb
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"
#include "frag_cvb.fh"
#include "formats_cvb.fh"
#include "WrkSpc.fh"

call rdioff_cvb(8,recinp,ioffs)
call rdr_cvb(work(ls(1)),nsyme*norb*norb,recinp,ioffs)
if ((ip(2) >= 1) .and. (.not. up2date_cvb('PRSYMELM'))) then
  do isyme=1,nsyme
    write(6,'(/,a,i4,3x,a)') ' Symmetry element no.',isyme,tags(isyme)
    ishft = norb*norb*(isyme-1)
    call mxprint_cvb(work(ishft+ls(1)),norb,norb,0)
  end do
  if (nsyme > 0) write(6,*) ' '
  call untouch_cvb('PRSYMELM')
end if

return

end subroutine mksymelm_cvb
