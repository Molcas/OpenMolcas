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

subroutine GetRest_t3(t1,t1_tmp,E2old)
! this file read 1) T1o
!                2) E1old,E2old,niter
! from RstFil file

use ChT3_global, only: LunAux, no, nv, printkey
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp), intent(out) :: t1(nv*no,2), t1_tmp(nv*no), E2old
integer(kind=iwp) :: dum
real(kind=wp) :: E1old

!open(unit=LunAux,File='RstFil',form='unformatted')
call MOLCAS_BinaryOpen_Vanilla(LunAux,'RstFil')
!mp write(u6,*) 'no, nv, length = ',no,nv,length
read(LunAux) t1(:,1)

call Map2_21_t3(t1,t1_tmp,nv,no)

t1(:,1) = t1_tmp
t1(:,2) = t1_tmp

read(LunAux) E1old,E2old,dum
#include "macros.fh"
unused_var(dum)

if (printkey > 1) write(u6,'(A,2(f15.12,1x))') 'Results from CCSD : E1, E2 ',E1old,E2old

close(LunAux)

return

end subroutine GetRest_t3
