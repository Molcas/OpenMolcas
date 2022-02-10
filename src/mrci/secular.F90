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

implicit real*8(A-H,O-Z)
intrinsic SQRT
dimension HMAT(NDIM,NDIM), SMAT(NDIM,NDIM)
dimension VEC(NDIM,NDIM), EVAL(NDIM), SCR(*)

THR2 = THR**2
! PUT NORMALIZED VECTORS INTO VEC:
call DCOPY_(N*NDIM,[0.0d00],0,VEC,1)
do I=1,N
  VEC(I,I) = 1.0d00/sqrt(SMAT(I,I))
end do
! GRAM-SCHMIDT ORTHONORMALIZING PROCEDURE:
NRON = 0
do I=1,N
  ! SMAT*(NORMALIZED VECTOR) INTO SCR:
  call DCOPY_(N,SMAT(1,I),1,SCR,1)
  call DSCAL_(N,VEC(I,I),SCR,1)
  ! PROJECT AWAY THE ALREADY ORTHONORMALIZED BASIS SET:
  do J=1,NRON
    MAXLEN = I-1-NRON+J
    SUM = 0.0d00
    do K=1,MAXLEN
      SUM = SUM+VEC(K,J)*SCR(K)
    end do
    do K=1,MAXLEN
      VEC(K,I) = VEC(K,I)-SUM*VEC(K,J)
    end do
  end do
  ! NORMALIZE AND MOVE INTO POSITION:
  SUM = 0.0d00
  do K=1,I
    SUM = SUM+VEC(K,I)*SCR(K)
  end do
  if (SUM < THR2) goto 60
  NRON = NRON+1
  SCALE = 1.0d00/sqrt(SUM)
  do K=1,I
    VEC(K,NRON) = SCALE*VEC(K,I)
  end do
60 continue
end do
do I=NRON+1,N
  call DCOPY_(N,[0.0d00],0,VEC(1,I),1)
end do
! TRANSFORM HAMILTONIAN INTO SCR:
IOFF1 = N*NRON
call DGEMM_('N','N',N,NRON,N,1.0d0,HMAT,NDIM,VEC,NDIM,0.0d0,SCR,N)
call DGEMM_('T','N',NRON,NRON,N,1.0d0,VEC,NDIM,SCR,N,0.0d0,SCR(IOFF1+1),NRON)
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
