!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 1986,1995, Bernd Artur Hess                            *
!               2005, Jesper Wisborg Krogh                             *
!***********************************************************************

subroutine DKRE1R(A,R,TT,V,G,RE1R,VEXTT,PVPT,N)
! CONSTRUCT RE1R FOR DK3

use Constants, only: Four
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: N
real(kind=wp), intent(in) :: A(N), R(N), TT(N), VEXTT(N*(N+1)/2), PVPT(N*(N+1)/2)
real(kind=wp), intent(out) :: V(N*(N+1)/2), G(N*(N+1)/2), RE1R(N,N)
integer(kind=iwp) :: I, IJ, J

IJ = 0
do I=1,N
  do J=1,I
    IJ = IJ+1
    V(IJ) = VEXTT(IJ)
    G(IJ) = PVPT(IJ)
  end do
end do

! MULTIPLY

IJ = 0
do I=1,N
  do J=1,I
    IJ = IJ+1
    V(IJ) = V(IJ)*A(I)*A(J)*R(I)*R(I)*R(J)*R(J)*TT(I)*TT(J)*Four
    RE1R(I,J) = V(IJ)
    RE1R(J,I) = RE1R(I,J)
  end do
end do

IJ = 0
do I=1,N
  do J=1,I
    IJ = IJ+1
    G(IJ) = G(IJ)*A(I)*A(J)*R(I)*R(J)
    RE1R(I,J) = RE1R(I,J)+G(IJ)
    RE1R(J,I) = RE1R(I,J)
  end do
end do

return

end subroutine DKRE1R
