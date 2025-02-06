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

subroutine GTJK_RASSCF(RJ,RK,NAC,IREOST)
! PURPOSE: GET ALL INTEGRALS COULOMB AND EXCHANGE INTEGRALS
!          WITH THE CHARGE DISTRIBUTION JK
!
! Input
! IREOST : Reorder array, symmetry => type (sic!)
!
! Modified by addition of IREOST, August 2003.

use wadr, only: TUVX
use Definitions, only: wp, iwp

implicit none
real(kind=wp) :: RJ(*), RK(*)
integer(kind=iwp) :: NAC, IREOST(*)
integer(kind=iwp) :: NT, NT_REO, NTT, NTU_REO, NTUJ, NTUK, NTUT, NU, NU_REO, NUT_REO

! FORM THE COULOMB (RJ) AND EXCHANGE (RK) INTEGRAL MATRICES FROM
! THE TWO-ELECTRON INTEGRAL LIST

NTUT = 0
do NT=1,NAC
  do NU=1,NT
    NT_REO = IREOST(NT)
    NU_REO = IREOST(NU)

    NTU_REO = NAC*(NT_REO-1)+NU_REO
    NUT_REO = NAC*(NU_REO-1)+NT_REO
    NTUT = NTUT+1
    NTUK = (NTUT**2+NTUT)/2
    RK(NTU_REO) = TUVX(NTUK)
    RK(NUT_REO) = TUVX(NTUK)

    NTT = (NT**2+NT)/2
    NTUJ = (NTT**2-NTT)/2+(NU**2+NU)/2
    RJ(NTU_REO) = TUVX(NTUJ)
    RJ(NUT_REO) = TUVX(NTUJ)
  end do
end do

end subroutine GTJK_RASSCF
