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

subroutine Even2r(idbg,N,V,G,E,A,R,TT,AUXF,AUXG,AUXH,W1W1)
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
integer(kind=iwp), intent(in) :: idbg, N
real(kind=wp), intent(inout) :: V(N*(N+1)/2), G(N*(N+1)/2)
real(kind=wp), intent(in) :: E(N), A(N), R(N), TT(N)
real(kind=wp), intent(out) :: AUXF(N,N), AUXG(N,N), AUXH(N,N), W1W1(N,N)
integer(kind=iwp) :: I, IE, IJ, J, M

#ifndef _DEBUGPRINT_
#include "macros.fh"
unused_var(idbg)
#endif

!ulf
#ifdef _DEBUGPRINT_
if (idbg > 0) then
  call PRMAT(IDBG,V,N,0,'V       ')
  call PRMAT(IDBG,G,N,0,'G       ')
  call PRMAT(IDBG,E,N,1,'E       ')
  call PRMAT(IDBG,A,N,1,'A       ')
  call PRMAT(IDBG,R,N,1,'R       ')
  call PRMAT(IDBG,TT,N,1,'TT      ')
end if
#endif
M = N
AUXH(:,:) = Zero
IJ = 0
do I=1,N
  do J=1,I
    IJ = IJ+1
    V(IJ) = V(IJ)/(E(I)+E(J))
    G(IJ) = G(IJ)/(E(I)+E(J))
  end do
end do
do J=1,N
  IJ = J*(J-1)/2+1
  do I=J,N
    IJ = IJ+I-1
    AUXF(I,J) = A(I)*R(I)*G(IJ)*A(J)*A(J)
    AUXG(I,J) = R(I)*V(IJ)*A(J)
  end do
end do
IJ = 0
do I=1,N
  do J=1,I
    IJ = IJ+1
    AUXF(J,I) = A(J)*R(J)*G(IJ)*A(I)*A(I)
    AUXG(J,I) = R(J)*V(IJ)*A(I)
  end do
end do

! ARQA ARQA

call CpLabr(AUXF,AUXG,N,N,N,M,M,AUXH,M,IE)
#ifdef _DEBUGPRINT_
if (IE /= 0) call SysHalt('relsew')
!ulf
if (idbg > 0) call prsq(idbg,'AUXH   1',auxh,n)
#endif
IJ = 0
do I=1,N
  do J=1,I
    IJ = IJ+1
    AUXG(J,I) = -Half/TT(J)*G(IJ)*A(I)*R(I)
  end do
end do
do J=1,N
  IJ = J*(J-1)/2+1
  do I=J,N
    IJ = IJ+I-1
    AUXG(I,J) = -Half/TT(I)*G(IJ)*A(J)*R(J)
  end do
end do

! ARQA AQRA

call CpLabr(AUXF,AUXG,N,N,N,M,M,AUXH,M,IE)
!ulf
#ifdef _DEBUGPRINT_
if (idbg > 0) call prsq(idbg,'AUXH   2',auxh,n)
#endif
IJ = 0
do I=1,N
  do J=1,I
    IJ = IJ+1
    AUXF(J,I) = A(J)*V(IJ)*A(I)*A(I)*R(I)
  end do
end do
do J=1,N
  IJ = J*(J-1)/2+1
  do I=J,N
    IJ = IJ+I-1
    AUXF(I,J) = A(I)*V(IJ)*A(J)*A(J)*R(J)
  end do
end do
IJ = 0
do I=1,N
  do J=1,I
    IJ = IJ+1
    AUXG(J,I) = -Two*TT(J)*R(J)*V(IJ)*A(I)
  end do
end do
do J=1,N
  IJ = J*(J-1)/2+1
  do I=J,N
    IJ = IJ+I-1
    AUXG(I,J) = -Two*TT(I)*R(I)*V(IJ)*A(J)
  end do
end do

! AQRA ARQA

call CpLabr(AUXF,AUXG,N,N,N,M,M,AUXH,M,IE)
!ulf
#ifdef _DEBUGPRINT_
if (idbg > 0) call prsq(idbg,'AUXH   3',auxh,n)
#endif
IJ = 0
do I=1,N
  do J=1,I
    IJ = IJ+1
    AUXG(J,I) = G(IJ)*A(I)*R(I)
  end do
end do
do J=1,N
  IJ = J*(J-1)/2+1
  do I=J,N
    IJ = IJ+I-1
    AUXG(I,J) = G(IJ)*A(J)*R(J)
  end do
end do

! AQRA AQRA

call CpLabr(AUXF,AUXG,N,N,N,M,M,AUXH,M,IE)

! KEEP W1*W1 FOR HIGHER-ORDER DK

W1W1(:,:) = AuxH

!ulf
#ifdef _DEBUGPRINT_
if (idbg > 0) call prsq(IDBG,'W*W     ',AUXH,N)
#endif

! 1/2 EW*W + 1/2 W*WE

do I=1,N
  do J=1,N
    AUXH(I,J) = Half*(AUXH(I,J)*E(I)+AUXH(I,J)*E(J))
  end do
end do
!ulf
#ifdef _DEBUGPRINT_
if (idbg > 0) call prsq(idbg,'AUXH SYM',auxh,n)
#endif

IJ = 0
do I=1,N
  do J=1,I
    IJ = IJ+1
    AUXF(J,I) = A(J)*R(J)*G(IJ)*A(I)*E(I)*A(I)
    AUXG(J,I) = R(J)*V(IJ)*A(I)
  end do
end do
do J=1,N
  IJ = J*(J-1)/2+1
  do I=J,N
    IJ = IJ+I-1
    AUXF(I,J) = A(I)*R(I)*G(IJ)*A(J)*E(J)*A(J)
    AUXG(I,J) = R(I)*V(IJ)*A(J)
  end do
end do
call CpLabr(AUXF,AUXG,N,N,N,M,M,AUXH,M,IE)
!ulf
#ifdef _DEBUGPRINT_
if (idbg > 0) call prsq(idbg,'AUXH   5',auxh,n)
#endif
IJ = 0
do I=1,N
  do J=1,I
    IJ = IJ+1
    AUXG(J,I) = -Half/TT(J)*G(IJ)*A(I)*R(I)
  end do
end do
do J=1,N
  IJ = J*(J-1)/2+1
  do I=J,N
    IJ = IJ+I-1
    AUXG(I,J) = -Half/TT(I)*G(IJ)*A(J)*R(J)
  end do
end do
call CpLabr(AUXF,AUXG,N,N,N,M,M,AUXH,M,IE)
!ulf
#ifdef _DEBUGPRINT_
if (idbg > 0) call prsq(idbg,'AUXH   6',auxh,n)
#endif
IJ = 0
do I=1,N
  do J=1,I
    IJ = IJ+1
    AUXF(J,I) = A(J)*V(IJ)*R(I)*A(I)*E(I)*A(I)
  end do
end do
do J=1,N
  IJ = J*(J-1)/2+1
  do I=J,N
    IJ = IJ+I-1
    AUXF(I,J) = A(I)*V(IJ)*R(J)*A(J)*E(J)*A(J)
  end do
end do
IJ = 0
do I=1,N
  do J=1,I
    IJ = IJ+1
    AUXG(J,I) = -Two*TT(J)*R(J)*V(IJ)*A(I)
  end do
end do
do J=1,N
  IJ = J*(J-1)/2+1
  do I=J,N
    IJ = IJ+I-1
    AUXG(I,J) = -Two*TT(I)*R(I)*V(IJ)*A(J)
  end do
end do
call CpLabr(AUXF,AUXG,N,N,N,M,M,AUXH,M,IE)
!ulf
#ifdef _DEBUGPRINT_
if (idbg > 0) call prsq(idbg,'AUXH   7',auxh,n)
#endif
IJ = 0
do I=1,N
  do J=1,I
    IJ = IJ+1
    AUXG(J,I) = G(IJ)*A(I)*R(I)
  end do
end do
do J=1,N
  IJ = J*(J-1)/2+1
  do I=J,N
    IJ = IJ+I-1
    AUXG(I,J) = G(IJ)*A(J)*R(J)
  end do
end do
call CpLabr(AUXF,AUXG,N,N,N,M,M,AUXH,M,IE)
!ulf
#ifdef _DEBUGPRINT_
if (idbg > 0) call prsq(idbg,'AUXH   8',auxh,n)
#endif

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

end subroutine Even2r
