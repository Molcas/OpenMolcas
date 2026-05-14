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
use MCLR_Data, only: Int2, NACOB
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(out) :: RJ(NACOB,NACOB), RK(NACOB,NACOB)
integer(kind=iwp) :: NT, NTUJ, NTUK, NU

! FORM THE COULOMB (RJ) AND EXCHANGE (RK) INTEGRAL MATRICES FROM
! THE TWO-ELECTRON INTEGRAL LIST

do NT=1,NACOB
  do NU=1,NT
    NTUK = iTri(nTri_Elem(nt),nTri_Elem(nu))
    RJ(NT,NU) = INT2(NTUK)
    RJ(NU,NT) = INT2(NTUK)

    NTUJ = nTri_Elem(iTri(nt,nu))
    RK(NT,NU) = INT2(NTUJ)
    RK(NU,NT) = INT2(NTUJ)
  end do
end do

end subroutine GTJK_MCLR
