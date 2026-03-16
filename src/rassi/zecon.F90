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

subroutine ZECON(NSTATE,N,UR,UI,AR,AI,ZEKL,IXYZ,ISTATE,ISS,JSS)

use definitions, only: iwp, wp

implicit none
integer(kind=iwp), intent(in) :: NSTATE, N, IXYZ, ISTATE, ISS, JSS
real(kind=wp), intent(in) :: UR(N,N), UI(N,N)
real(kind=wp), intent(in) :: AR(N,N), AI(N,N)
complex(kind=wp), intent(inout) :: ZEKL(2,2,3,NSTATE)
real(kind=wp) TMPR1, TMPR2, TMPI1, TMPI2

TMPR1 = AR(ISS,JSS)*UR(JSS,1)-AI(ISS,JSS)*UI(JSS,1)
TMPR2 = AR(ISS,JSS)*UR(JSS,2)-AI(ISS,JSS)*UI(JSS,2)
TMPI1 = AI(ISS,JSS)*UR(JSS,1)+AR(ISS,JSS)*UI(JSS,1)
TMPI2 = AI(ISS,JSS)*UR(JSS,2)+AR(ISS,JSS)*UI(JSS,2)

ZEKL(1,1,IXYZ,ISTATE) = ZEKL(1,1,IXYZ,ISTATE)+cmplx(UR(ISS,1)*TMPR1+UI(ISS,1)*TMPI1,UR(ISS,1)*TMPI1-UI(ISS,1)*TMPR1,kind=wp)
ZEKL(1,2,IXYZ,ISTATE) = ZEKL(1,2,IXYZ,ISTATE)+cmplx(UR(ISS,1)*TMPR2+UI(ISS,1)*TMPI2,UR(ISS,1)*TMPI2-UI(ISS,1)*TMPR2,kind=wp)
ZEKL(2,1,IXYZ,ISTATE) = ZEKL(2,1,IXYZ,ISTATE)+cmplx(UR(ISS,2)*TMPR1+UI(ISS,2)*TMPI1,UR(ISS,2)*TMPI1-UI(ISS,2)*TMPR1,kind=wp)
ZEKL(2,2,IXYZ,ISTATE) = ZEKL(2,2,IXYZ,ISTATE)+cmplx(UR(ISS,2)*TMPR2+UI(ISS,2)*TMPI2,UR(ISS,2)*TMPI2-UI(ISS,2)*TMPR2,kind=wp)

end subroutine ZECON
