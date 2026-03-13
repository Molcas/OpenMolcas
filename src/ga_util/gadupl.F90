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
! Copyright (C) 1995, Martin Schuetz                                   *
!               1998, Roland Lindh                                     *
!               2000-2015, Steven Vancoillie                           *
!***********************************************************************

! duplicate & copy a global array...
! iGA1,iGA2:       GA handles...
subroutine GADupl(iGA1,iGA2)

#ifdef _MOLCAS_MPP_
use Para_Info, only: Is_Real_Par
use Definitions, only: u6
#endif
use Definitions, only: iwp

implicit none
#ifdef _MOLCAS_MPP_
integer(kind=iwp), intent(in) :: iGA1
integer(kind=iwp), intent(out) :: iGA2
logical(kind=iwp) :: ok
character(len=6) :: gaLbl2
character(len=5) :: gaLbl
#include "global.fh"

if (.not. Is_Real_Par()) return
if (iGA1 >= 0) return
call ga_inquire_name(iGA1,gaLbl)
write(gaLbl2,'(A,I1)') gaLbl,2

ok = ga_duplicate(iGA1,iGA2,gaLbl2)
if (.not. ok) then
  write(u6,*) 'GADupl: ga_duplicate not OK!'
  call ga_error('GADupl',42)
  call Abend()
end if

call ga_copy(iGA1,iGA2)
#else
integer(kind=iwp), intent(in) :: iGA1, iGA2
#include "macros.fh"
unused_var(iGA1)
unused_var(iGA2)
#endif

end subroutine GADupl
