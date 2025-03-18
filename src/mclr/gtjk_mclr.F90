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

use Arrays, only: Int2
use MCLR_Data, only: NACOB

implicit none
real*8 RJ(NACOB,NACOB), RK(NACOB,nACOB)
integer NT, NU, NTUK, NTUJ
integer i, j, itri
! Statement function
itri(i,j) = max(i,j)*(max(i,j)-1)/2+min(i,j)

! FORM THE COULOMB (RJ) AND EXCHANGE (RK) INTEGRAL MATRICES FROM
! THE TWO-ELECTRON INTEGRAL LIST

do NT=1,NACOB
  do NU=1,NT
    NTUK = itri(itri(nt,nt),itri(nu,nu))
    RJ(NT,NU) = INT2(NTUK)
    RJ(NU,NT) = INT2(NTUK)

    NTUJ = itri(itri(nt,nu),itri(nu,nt))
    RK(NT,NU) = INT2(NTUJ)
    RK(NU,NT) = INT2(NTUJ)
  end do
end do

end subroutine GTJK_MCLR
