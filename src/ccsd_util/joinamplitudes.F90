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

subroutine joinamplitudes(wrk,wrksize)
! this routine joins contributions to amplitudes (Tn) from all
! processors. Sum will be placed back to Tn on all machines
!
! N.B. FREE ZONE in purpose of this routine is the space of
! free working files V1-V4, .....
! (free zone is not used in GA)
! Since T2 amplitudes are no more then oovv, at most V1-V3 space
! will be damaged. (actually V1 and V2 only)

use Para_Info, only: nProcs
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: wrksize
real(kind=wp) :: wrk(wrksize)
#include "ccsd2.fh"
integer(kind=iwp) :: ii, length

if (nProcs == 1) return

!1 join t13,t14

!1.1 calc overall length of t13 and t14 (they must be one after the other)
ii = mapdt14(0,5)
length = mapdt14(ii,1)+mapdt14(ii,2)-posst130

!1.2 vanish required part in free zone
!call mv0zero(length,length,wrk(possv10))

!1.3 allreduce t13 and t14 into free zone
!call MPI_ALLREDUCE(wrk(posst130),wrk(possv10),length,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,rc)

!1.4 put joined t13,t14 back to t13,t14 place from free zone
!do i=0,length-1
!  wrk(posst130+i) = wrk(possv10+i)
!end do

!1.om allreduce t13,t14 together
call gadgop(wrk(posst130),length,'+')

!2 join t21,t22,t23

!2.1 calc overall length of t21-t23 (they must be one after the other)
ii = mapdt23(0,5)
length = mapdt23(ii,1)+mapdt23(ii,2)-posst210

!2.2 vanish required part in free zone
!call mv0zero(length,length,wrk(possv10))

!2.3 allreduce t13 and t14 into free zone
!call MPI_ALLREDUCE (wrk(posst210),wrk(possv10),length,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,rc)

!2.4 put joined t21-t23 back to t21-t23 place from free zone
!do i=0,length-1
!  wrk(posst210+i) = wrk(possv10+i)
!end do

!2.om allreduce t21-t23 together
call gadgop(wrk(posst210),length,'+')

return

end subroutine joinamplitudes
