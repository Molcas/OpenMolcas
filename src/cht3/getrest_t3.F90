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

use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp) :: t1(*), t1_tmp(*), E2old
#include "cht3_ccsd1.fh"
#include "ccsd_t3compat.fh"
integer(kind=iwp) :: i, length
real(kind=wp) :: E1old

!open(unit=LunAux,File='RstFil',form='unformatted')
call MOLCAS_BinaryOpen_Vanilla(LunAux,'RstFil')
length = nv*no
!mp write(u6,*) 'no, nv, length = ',no,nv,length
read(LunAux) t1(1:length)

call Map2_21_t3(t1,t1_tmp,nv,no)

do i=1,length
  t1(i+length) = t1_tmp(i)
  t1(i) = t1_tmp(i)
end do

read(LunAux) E1old,E2old,i

if (printkey > 1) write(u6,'(A,2(f15.12,1x))') 'Results from CCSD : E1, E2 ',E1old,E2old

close(LunAux)

return

end subroutine GetRest_t3
