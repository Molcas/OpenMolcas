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

use ccsd_global, only: t13, t14, t21, t23
use Para_Info, only: nProcs
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: wrksize
real(kind=wp), intent(inout) :: wrk(wrksize)
integer(kind=iwp) :: ii, length

if (nProcs == 1) return

!1 join t13,t14

!1.1 calc overall length of t13 and t14 (they must be one after the other)
ii = t14%d(0,5)
length = t14%d(ii,1)+t14%d(ii,2)-t13%pos0

!1.2 vanish required part in free zone
!wrk(v1%pos0:v1%pos0+length-1) = Zero

!1.3 allreduce t13 and t14 into free zone
!call MPI_ALLREDUCE(wrk(t13%pos0),wrk(v1%pos0),length,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,rc)

!1.4 put joined t13,t14 back to t13,t14 place from free zone
!do i=0,length-1
!  wrk(t13%pos0+i) = wrk(v1%pos0+i)
!end do

!1.om allreduce t13,t14 together
call gadgop(wrk(t13%pos0),length,'+')

!2 join t21,t22,t23

!2.1 calc overall length of t21-t23 (they must be one after the other)
ii = t23%d(0,5)
length = t23%d(ii,1)+t23%d(ii,2)-t21%pos0

!2.2 vanish required part in free zone
!wrk(v1%pos0:v1%pos0+length-1) = Zero

!2.3 allreduce t13 and t14 into free zone
!call MPI_ALLREDUCE (wrk(t21%pos0),wrk(v1%pos0),length,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,rc)

!2.4 put joined t21-t23 back to t21-t23 place from free zone
!do i=0,length-1
!  wrk(t21%pos0+i) = wrk(v1%pos0+i)
!end do

!2.om allreduce t21-t23 together
call gadgop(wrk(t21%pos0),length,'+')

return

end subroutine joinamplitudes
