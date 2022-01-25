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

implicit real*8(A-H,O-Z)
dimension V(N*(N+1)/2), G(N*(N+1)/2), E(N), R(N), A(N), TT(N), AUXF(N,N), AUXG(N,N), AUXH(N,N)
dimension EVEN1(N,N)
dimension RE1R(N,N)
dimension VEXTT(N*(N+1)/2), PVPT(N*(N+1)/2)
dimension AUXI(N,N)
dimension W1W1(N,N)

! CONSTRUCT RE1R

call DKRE1R(A,R,E,TT,V,G,RE1R,VEXTT,PVPT,N)

M = N
IJ = 0
do I=1,N
  do J=1,I
    IJ = IJ+1
    AUXH(I,J) = 0.d0
    AUXH(J,I) = 0.d0
    V(IJ) = VEXTT(IJ)/(E(I)+E(J))
    G(IJ) = PVPT(IJ)/(E(I)+E(J))
  end do
end do

! W1*W1*E1

! 1/2 E1*W1*W1 + 1/2 W1*W1*E1

do I=1,N
  do J=1,N
    AUXF(I,J) = 0.5d0*EVEN1(I,J)
  end do
end do
call CpLabr(W1W1,AUXF,N,N,N,M,M,AUXH,M,IE)
call CpLabr(AUXF,W1W1,N,N,N,M,M,AUXH,M,IE)

! W1*E1*W1 TERM

IJ = 0
do I=1,N
  do J=1,I
    IJ = IJ+1
    AUXF(I,J) = A(I)*R(I)*G(IJ)*A(J)/R(J)/TT(J)*0.5d0
    AUXF(J,I) = A(J)*R(J)*G(IJ)*A(I)/R(I)/TT(I)*0.5d0
    AUXG(I,J) = -A(I)*V(IJ)*A(J)
    AUXG(J,I) = -A(J)*V(IJ)*A(I)
  end do
end do
call CZERO2(AUXI,N,N,N)
call CpLabr(AUXF,RE1R,N,N,N,M,M,AUXI,M,IE)
call CpLabr(AUXI,AUXG,N,N,N,M,M,AUXH,M,IE)

IJ = 0
do I=1,N
  do J=1,I
    IJ = IJ+1
    AUXF(I,J) = A(I)*V(IJ)*A(J)
    AUXF(J,I) = A(J)*V(IJ)*A(I)
    AUXG(I,J) = -A(I)/R(I)*G(IJ)*A(J)*R(J)/TT(I)*0.5d0
    AUXG(J,I) = -A(J)/R(J)*G(IJ)*A(I)*R(I)/TT(J)*0.5d0
  end do
end do
call CZERO2(AUXI,N,N,N)
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
call CZERO2(AUXI,N,N,N)
call CpLabr(AUXF,RE1R,N,N,N,M,M,AUXI,M,IE)
call CpLabr(AUXI,AUXG,N,N,N,M,M,AUXH,M,IE)

IJ = 0
do I=1,N
  do J=1,I
    IJ = IJ+1
    AUXF(I,J) = A(I)*R(I)*G(IJ)*A(J)/R(J)/TT(J)*0.5d0
    AUXF(J,I) = A(J)*R(J)*G(IJ)*A(I)/R(I)/TT(I)*0.5d0
    AUXG(I,J) = A(I)/R(I)*G(IJ)*A(J)*R(J)/TT(I)*0.5d0
    AUXG(J,I) = A(J)/R(J)*G(IJ)*A(I)*R(I)/TT(J)*0.5d0
  end do
end do
call CZERO2(AUXI,N,N,N)
call CpLabr(AUXF,RE1R,N,N,N,M,M,AUXI,M,IE)
call CpLabr(AUXI,AUXG,N,N,N,M,M,AUXH,M,IE)

! SYMMETRISIEREN

IJ = 0
do I=1,N
  do J=1,I
    IJ = IJ+1
    G(IJ) = 0.5d0*(AUXH(I,J)+AUXH(J,I))
  end do
end do

return

end subroutine EVEN3
