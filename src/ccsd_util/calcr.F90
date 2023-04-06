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

subroutine calcr(wrk,wrksize,lune)
! this routine calcs difference vector Tn = Tn-E
! Tn=(T12,T22,T23,T13,T14)
!
! lune - lun of file, where E is stored (I)

use ccsd_global, only: t13, t14, t21, t22, t23, v1
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: wrksize
integer(kind=iwp), intent(_IN_) :: lune
real(kind=wp), intent(inout) :: wrk(wrksize)
integer(kind=iwp) :: rc

!1 rewind lune
call filemanager(2,lune,rc)

! T2aaaa
call getmediate(wrk,wrksize,lune,v1,rc)
call calcrh1(wrk,wrksize,t21,v1)
! T2bbbb
call getmediate(wrk,wrksize,lune,v1,rc)
call calcrh1(wrk,wrksize,t22,v1)
! T2abab
call getmediate(wrk,wrksize,lune,v1,rc)
call calcrh1(wrk,wrksize,t23,v1)
! T1aa
call getmediate(wrk,wrksize,lune,v1,rc)
call calcrh1(wrk,wrksize,t13,v1)
! T1bb
call getmediate(wrk,wrksize,lune,v1,rc)
call calcrh1(wrk,wrksize,t14,v1)

return

end subroutine calcr
