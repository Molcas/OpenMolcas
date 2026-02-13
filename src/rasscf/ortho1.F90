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

subroutine ORTHO1(S,U,V,N,M)
! Purpose: Orthonormalize N times N vector set U.
!
! Called from ORTHO.
!
! Subroutine calls: ORTHO2.
!
! ********** IBM-3090 Release 88 10 10 **********

use Constants, only: One
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(in) :: S(*)
real(kind=wp), intent(inout) :: U(*), V(*)
integer(kind=iwp), intent(in) :: N, M
real(kind=wp) XNORM, OVL, OVL1
real(kind=wp), parameter :: THR = 0.2_wp
real(kind=wp), external :: DDot_
integer(kind=iwp) :: IBASE, I, JBASE, J

IBASE = 1

do I=1,M

  do
    call ORTHO2(S,U(IBASE),V(IBASE),N)
    ! Normalize U and calculate V=S*U.
    XNORM = ONE
    JBASE = 1
    do J=1,I-1
      OVL = DDOT_(N,U(IBASE),1,V(JBASE),1)
      OVL1 = -OVL
      call DAXPY_(N,OVL1,U(JBASE),1,U(IBASE),1)
      XNORM = XNORM-OVL**2
      if (XNORM >= THR) exit
      JBASE = JBASE+N
    end do
    call ORTHO2(S,U(IBASE),V(IBASE),N)
    IBASE = IBASE+N
  end do

end do

end subroutine ORTHO1
