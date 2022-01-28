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
! Copyright (C) 1986, Bernd Artur Hess                                 *
!***********************************************************************

subroutine EVEN2(N,V,G,E,A,R,TT,AUXF,AUXG,AUXH,W1W1)
!     EVEN2 - BERND HESS - V 1.0 - 5.2.86
!     CALCULATE EVEN2 OPERATORS
!
!     N       DIMENSION OF MATRICES
!     V       POTENTIAL MATRIX
!     G       MATRIX OF PVP OPERATOR. WILL CONTAIN EVEN2 OPERATORS
!             ON OUTPUT
!     E       RELATIVISTIC ENERGY (DIAGONAL)
!     A       A-FACTORS (DIAGONAL)
!     R       R-FACTORS (DIAGONAL)
!     TT      NONREL. KINETIC ENERGY (DIAGONAL)
!     AUXF,AUXG,AUXH  SCRATCH ARAYS

use Constants, only: Zero, Two, Half
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: N
real(kind=wp), intent(inout) :: V(N*(N+1)/2), G(N*(N+1)/2)
real(kind=wp), intent(in) :: E(N), A(N), R(N), TT(N)
real(kind=wp), intent(out) :: AUXF(N,N), AUXG(N,N), AUXH(N,N), W1W1(N,N)
integer(kind=iwp) :: I, IE, IJ, J, M

M = N
IJ = 0
do I=1,N
  do J=1,I
    IJ = IJ+1
    AUXH(I,J) = Zero
    AUXH(J,I) = Zero
    V(IJ) = V(IJ)/(E(I)+E(J))
    G(IJ) = G(IJ)/(E(I)+E(J))
    AUXF(I,J) = A(I)*R(I)*G(IJ)*A(J)*A(J)
    AUXF(J,I) = A(J)*R(J)*G(IJ)*A(I)*A(I)
    AUXG(I,J) = R(I)*V(IJ)*A(J)
    AUXG(J,I) = R(J)*V(IJ)*A(I)
  end do
end do

! ARQA ARQA

call CPLAB(AUXF,AUXG,N,N,N,M,M,AUXH,M,IE)
if (IE /= 0) call SysHalt('relint')
!call PRSQ('AUXH   1',AUXH,N)
IJ = 0
do I=1,N
  do J=1,I
    IJ = IJ+1
    AUXG(I,J) = -Half/TT(I)*G(IJ)*A(J)*R(J)
    AUXG(J,I) = -Half/TT(J)*G(IJ)*A(I)*R(I)
  end do
end do

! ARQA AQRA

call CPLAB(AUXF,AUXG,N,N,N,M,M,AUXH,M,IE)
!call PRSQ('AUXH   2',AUXH,N)
IJ = 0
do I=1,N
  do J=1,I
    IJ = IJ+1
    AUXF(I,J) = A(I)*V(IJ)*A(J)*A(J)*R(J)
    AUXF(J,I) = A(J)*V(IJ)*A(I)*A(I)*R(I)
    AUXG(I,J) = -Two*TT(I)*R(I)*V(IJ)*A(J)
    AUXG(J,I) = -Two*TT(J)*R(J)*V(IJ)*A(I)
  end do
end do

! AQRA ARQA

call CPLAB(AUXF,AUXG,N,N,N,M,M,AUXH,M,IE)
!call PRSQ('AUXH   3',AUXH,N)
IJ = 0
do I=1,N
  do J=1,I
    IJ = IJ+1
    AUXG(I,J) = G(IJ)*A(J)*R(J)
    AUXG(J,I) = G(IJ)*A(I)*R(I)
  end do
end do

! AQRA AQRA

call CPLAB(AUXF,AUXG,N,N,N,M,M,AUXH,M,IE)

! KEEP W1*W1 FOR HIGHER-ORDER DK

do I=1,N
  do J=1,N
    W1W1(I,J) = AUXH(I,J)
  end do
end do

! 1/2 EW*W + 1/2 W*WE

do I=1,N
  do J=1,N
    AUXH(I,J) = Half*(AUXH(I,J)*E(I)+AUXH(I,J)*E(J))
  end do
end do

IJ = 0
do I=1,N
  do J=1,I
    IJ = IJ+1
    AUXF(I,J) = A(I)*R(I)*G(IJ)*A(J)*E(J)*A(J)
    AUXF(J,I) = A(J)*R(J)*G(IJ)*A(I)*E(I)*A(I)
    AUXG(I,J) = R(I)*V(IJ)*A(J)
    AUXG(J,I) = R(J)*V(IJ)*A(I)
  end do
end do
call CPLAB(AUXF,AUXG,N,N,N,M,M,AUXH,M,IE)
IJ = 0
do I=1,N
  do J=1,I
    IJ = IJ+1
    AUXG(I,J) = -Half/TT(I)*G(IJ)*A(J)*R(J)
    AUXG(J,I) = -Half/TT(J)*G(IJ)*A(I)*R(I)
  end do
end do
call CPLAB(AUXF,AUXG,N,N,N,M,M,AUXH,M,IE)
!call PRSQ('AUXH   6',AUXH,N)
IJ = 0
do I=1,N
  do J=1,I
    IJ = IJ+1
    AUXF(I,J) = A(I)*V(IJ)*R(J)*A(J)*E(J)*A(J)
    AUXF(J,I) = A(J)*V(IJ)*R(I)*A(I)*E(I)*A(I)
    AUXG(I,J) = -Two*TT(I)*R(I)*V(IJ)*A(J)
    AUXG(J,I) = -Two*TT(J)*R(J)*V(IJ)*A(I)
  end do
end do
call CPLAB(AUXF,AUXG,N,N,N,M,M,AUXH,M,IE)
IJ = 0
do I=1,N
  do J=1,I
    IJ = IJ+1
    AUXG(I,J) = G(IJ)*A(J)*R(J)
    AUXG(J,I) = G(IJ)*A(I)*R(I)
  end do
end do
call CPLAB(AUXF,AUXG,N,N,N,M,M,AUXH,M,IE)
!call PRSQ('AUXH   8',AUXH,N)

! SYMMETRISIEREN

IJ = 0
do I=1,N
  do J=1,I
    IJ = IJ+1
    G(IJ) = -Half*(AUXH(I,J)+AUXH(J,I))
  end do
end do
!call PRM('OUTPUT  ',G,N)

return

end subroutine EVEN2
