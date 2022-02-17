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

subroutine SECULAR(NDIM,N,NRON,HMAT,SMAT,VEC,EVAL,SCR,THR)

use Constants, only: Zero, One
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: NDIM, N
integer(kind=iwp), intent(out) :: NRON
real(kind=wp), intent(in) :: HMAT(NDIM,NDIM), SMAT(NDIM,NDIM), THR
real(kind=wp), intent(out) :: VEC(NDIM,NDIM), EVAL(NDIM)
real(kind=wp), intent(_OUT_) :: SCR(*)
integer(kind=iwp) :: I, IFROM, II, IOFF1, ITO, J, K, MAXLEN
real(kind=wp) :: RSUM, SCL, THR2

THR2 = THR**2
! PUT NORMALIZED VECTORS INTO VEC:
VEC(:,1:N) = Zero
do I=1,N
  VEC(I,I) = One/sqrt(SMAT(I,I))
end do
! GRAM-SCHMIDT ORTHONORMALIZING PROCEDURE:
NRON = 0
do I=1,N
  ! SMAT*(NORMALIZED VECTOR) INTO SCR:
  SCR(1:N) = VEC(I,I)*SMAT(1:N,I)
  ! PROJECT AWAY THE ALREADY ORTHONORMALIZED BASIS SET:
  do J=1,NRON
    MAXLEN = I-1-NRON+J
    RSUM = Zero
    do K=1,MAXLEN
      RSUM = RSUM+VEC(K,J)*SCR(K)
    end do
    do K=1,MAXLEN
      VEC(K,I) = VEC(K,I)-RSUM*VEC(K,J)
    end do
  end do
  ! NORMALIZE AND MOVE INTO POSITION:
  RSUM = Zero
  do K=1,I
    RSUM = RSUM+VEC(K,I)*SCR(K)
  end do
  if (RSUM >= THR2) then
    NRON = NRON+1
    SCL = One/sqrt(RSUM)
    do K=1,I
      VEC(K,NRON) = SCL*VEC(K,I)
    end do
  end if
end do
do I=NRON+1,N
  VEC(1:N,I) = Zero
end do
! TRANSFORM HAMILTONIAN INTO SCR:
IOFF1 = N*NRON
call DGEMM_('N','N',N,NRON,N,One,HMAT,NDIM,VEC,NDIM,Zero,SCR,N)
call DGEMM_('T','N',NRON,NRON,N,One,VEC,NDIM,SCR,N,Zero,SCR(IOFF1+1),NRON)
! COPY TRANSFORMED HMAT INTO TRIANGULAR STORAGE IN SCR:
IFROM = IOFF1+1
ITO = 1
do I=1,NRON
  call DCOPY_(I,SCR(IFROM),1,SCR(ITO),1)
  ITO = ITO+I
  IFROM = IFROM+NRON
end do
! DIAGONALIZE:
call Jacob(SCR,VEC,NRON,NDIM)
! COPY EIGENVALUES INTO EVAL:
II = 0
do I=1,NRON
  II = II+I
  EVAL(I) = SCR(II)
end do

return

end subroutine SECULAR
