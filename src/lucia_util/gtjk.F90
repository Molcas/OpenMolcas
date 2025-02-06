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

subroutine GTJK(RJ,RK,NTOOB,IREOST)
! Interface routine for obtaining Coulomb (RJ) and
! Exchange integrals (RK)
!
! Ordering of integrals is in the internal order

use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: NTOOB, IREOST(*)
real(kind=wp) :: RJ(NTOOB,NTOOB), RK(NTOOB,NTOOB)
integer(kind=iwp) :: NTEST

call GTJK_RASSCF(RJ,RK,NTOOB,IREOST)

NTEST = 0
if (NTEST /= 0) then
  write(u6,*) ' RJ and RK from GTJK'
  call WRTMAT(RJ,NTOOB,NTOOB,NTOOB,NTOOB)
  call WRTMAT(RK,NTOOB,NTOOB,NTOOB,NTOOB)
end if

end subroutine GTJK
