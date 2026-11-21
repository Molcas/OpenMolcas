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

subroutine GTJK_MCLR(RJ,RK)
! PURPOSE: GET ALL INTEGRALS COULOMB AND EXCHANGE INTEGRALS
!          WITH THE CHARGE DISTRIBUTION JK

use Index_Functions, only: iTri, nTri_Elem
use MCLR_Data, only: Int2, IREOTS, NACOB
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(out) :: RJ(NACOB,NACOB), RK(NACOB,NACOB)
integer(kind=iwp) :: iT, iU, NT, NTUJ, NTUK, NU

! FORM THE COULOMB (RJ) AND EXCHANGE (RK) INTEGRAL MATRICES FROM
! THE TWO-ELECTRON INTEGRAL LIST

! RJ, RK in Lucia's (type) order of the orbitals, in which CIDIA4
! addresses them; INT2 is addressed with symmetry-ordered indices
do iT=1,NACOB
  NT = IREOTS(iT)
  do iU=1,iT
    NU = IREOTS(iU)
    NTUK = iTri(nTri_Elem(nt),nTri_Elem(nu))
    RJ(iT,iU) = INT2(NTUK)
    RJ(iU,iT) = INT2(NTUK)

    NTUJ = nTri_Elem(iTri(nt,nu))
    RK(iT,iU) = INT2(NTUJ)
    RK(iU,iT) = INT2(NTUJ)
  end do
end do

end subroutine GTJK_MCLR
