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

!#define _DEBUGPRINT_
subroutine GTJK(RJ,RK,NTOOB,IREOST)
! Interface routine for obtaining Coulomb (RJ) and
! Exchange integrals (RK)
!
! Ordering of integrals is in the internal order

use Index_Functions, only: nTri_Elem
use wadr, only: TUVX
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: NTOOB, IREOST(NTOOB)
real(kind=wp), intent(inout) :: RJ(NTOOB,NTOOB), RK(NTOOB,NTOOB)
integer(kind=iwp) :: NT, NT_REO, NTT, NTUJ, NTUK, NTUT, NU, NU_REO

NTUT = 0
do NT=1,NTOOB
  do NU=1,NT
    NT_REO = IREOST(NT)
    NU_REO = IREOST(NU)

    NTUT = NTUT+1
    NTUK = nTri_Elem(NTUT)
    RK(NT_REO,NU_REO) = TUVX(NTUK)
    RK(NU_REO,NT_REO) = TUVX(NTUK)

    NTT = nTri_Elem(NT)
    NTUJ = nTri_Elem(NTT-1)+nTri_Elem(NU)
    RJ(NT_REO,NU_REO) = TUVX(NTUJ)
    RJ(NU_REO,NT_REO) = TUVX(NTUJ)
  end do
end do

#ifdef _DEBUGPRINT_
write(u6,*) ' RJ and RK from GTJK'
call WRTMAT(RJ,NTOOB,NTOOB,NTOOB,NTOOB)
call WRTMAT(RK,NTOOB,NTOOB,NTOOB,NTOOB)
#endif

end subroutine GTJK
