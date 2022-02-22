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

subroutine EVEN3(N,V,G,E,A,R,TT,AUXF,AUXG,AUXH,EVEN1,VEXTT,PVPT,RE1R,W1W1,AUXI)

use Constants, only: Zero, Half
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: N
real(kind=wp), intent(out) :: V(N*(N+1)/2), G(N*(N+1)/2), AUXF(N,N), AUXG(N,N), AUXH(N,N), RE1R(N,N), AUXI(N,N)
real(kind=wp), intent(in) :: E(N), A(N), R(N), TT(N), EVEN1(N,N), VEXTT(N*(N+1)/2), PVPT(N*(N+1)/2), W1W1(N,N)
integer(kind=iwp) :: I, IE, IJ, J, M

! CONSTRUCT RE1R

call DKRE1R(A,R,TT,V,G,RE1R,VEXTT,PVPT,N)

M = N
IJ = 0
do I=1,N
  do J=1,I
    IJ = IJ+1
    AUXH(I,J) = Zero
    AUXH(J,I) = Zero
    V(IJ) = VEXTT(IJ)/(E(I)+E(J))
    G(IJ) = PVPT(IJ)/(E(I)+E(J))
  end do
end do

! W1*W1*E1

! 1/2 E1*W1*W1 + 1/2 W1*W1*E1

do I=1,N
  do J=1,N
    AUXF(I,J) = Half*EVEN1(I,J)
  end do
end do
call CpLabr(W1W1,AUXF,N,N,N,M,M,AUXH,M,IE)
call CpLabr(AUXF,W1W1,N,N,N,M,M,AUXH,M,IE)

! W1*E1*W1 TERM

IJ = 0
do I=1,N
  do J=1,I
    IJ = IJ+1
    AUXF(I,J) = A(I)*R(I)*G(IJ)*A(J)/R(J)/TT(J)*Half
    AUXF(J,I) = A(J)*R(J)*G(IJ)*A(I)/R(I)/TT(I)*Half
    AUXG(I,J) = -A(I)*V(IJ)*A(J)
    AUXG(J,I) = -A(J)*V(IJ)*A(I)
  end do
end do
AUXI(:,:) = Zero
call CpLabr(AUXF,RE1R,N,N,N,M,M,AUXI,M,IE)
call CpLabr(AUXI,AUXG,N,N,N,M,M,AUXH,M,IE)

IJ = 0
do I=1,N
  do J=1,I
    IJ = IJ+1
    AUXF(I,J) = A(I)*V(IJ)*A(J)
    AUXF(J,I) = A(J)*V(IJ)*A(I)
    AUXG(I,J) = -A(I)/R(I)*G(IJ)*A(J)*R(J)/TT(I)*Half
    AUXG(J,I) = -A(J)/R(J)*G(IJ)*A(I)*R(I)/TT(J)*Half
  end do
end do
AUXI(:,:) = Zero
call CpLabr(AUXF,RE1R,N,N,N,M,M,AUXI,M,IE)
call CpLabr(AUXI,AUXG,N,N,N,M,M,AUXH,M,IE)

IJ = 0
do I=1,N
  do J=1,I
    IJ = IJ+1
    AUXF(I,J) = A(I)*V(IJ)*A(J)
    AUXF(J,I) = A(J)*V(IJ)*A(I)
    AUXG(I,J) = A(I)*V(IJ)*A(J)
    AUXG(J,I) = A(J)*V(IJ)*A(I)
  end do
end do
AUXI(:,:) = Zero
call CpLabr(AUXF,RE1R,N,N,N,M,M,AUXI,M,IE)
call CpLabr(AUXI,AUXG,N,N,N,M,M,AUXH,M,IE)

IJ = 0
do I=1,N
  do J=1,I
    IJ = IJ+1
    AUXF(I,J) = A(I)*R(I)*G(IJ)*A(J)/R(J)/TT(J)*Half
    AUXF(J,I) = A(J)*R(J)*G(IJ)*A(I)/R(I)/TT(I)*Half
    AUXG(I,J) = A(I)/R(I)*G(IJ)*A(J)*R(J)/TT(I)*Half
    AUXG(J,I) = A(J)/R(J)*G(IJ)*A(I)*R(I)/TT(J)*Half
  end do
end do
AUXI(:,:) = Zero
call CpLabr(AUXF,RE1R,N,N,N,M,M,AUXI,M,IE)
call CpLabr(AUXI,AUXG,N,N,N,M,M,AUXH,M,IE)

! SYMMETRISIEREN

IJ = 0
do I=1,N
  do J=1,I
    IJ = IJ+1
    G(IJ) = Half*(AUXH(I,J)+AUXH(J,I))
  end do
end do

return

end subroutine EVEN3
