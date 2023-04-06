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

subroutine saverest1(wrk,wrksize,lunrst)
! this routine saves restart informations:
! t13,t14,t21,t22,t23

use ccsd_global, only: keyrst, t13, t14, t21, t22, t23
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: wrksize
real(kind=wp), intent(_IN_) :: wrk(wrksize)
integer(kind=iwp), intent(_IN_) :: lunrst
integer(kind=iwp) :: rc

!0 return if need
if (keyrst == 0) return

!1 rewind tape
call filemanager(2,lunrst,rc)

!2 write T1aa
call wrtmediate(wrk,wrksize,lunrst,t13,rc)

!3 write T1bb
call wrtmediate(wrk,wrksize,lunrst,t14,rc)

!4 write T2aaaa
call wrtmediate(wrk,wrksize,lunrst,t21,rc)

!5 write T2bbbb
call wrtmediate(wrk,wrksize,lunrst,t22,rc)

!6 write T2abab
call wrtmediate(wrk,wrksize,lunrst,t23,rc)

return

end subroutine saverest1
