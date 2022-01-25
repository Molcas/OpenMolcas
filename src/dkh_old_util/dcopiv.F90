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
! Copyright (C) 1988, H. Rieger                                        *
!***********************************************************************
!**********************************************************************
!
!        PROGRAMMBIBLIOTHEK RHRZ BONN        23/11/88       DCOPIV
!                                            FORTRAN 77     IBM 3081
!
!
! NAME:    DCOPIV
!
! PURPOSE:
!
! This program calculates the solution of the system of linear
! equations  A*X = B  for several right-hand sides. The method used is
! complete pivoting which is generally considered a stable method,
! even if the matrix A is ill-conditioned.
! A is decomposed into the product of a lower triangular matrx L
! with main diagonal elements equal to 1 and an upper triangular
! matrix U such that:  L*U = A. While L is not saved, U is over-
! written on A.
! If the matrix A is nonsingular, the solutions are stored in B.
!
! USAGE:   CALL DCOPIV(A,B,N,M,D,EPS,DET,EX,CTR,S)
!
! PARAMETERS:
!
! N,M:     INTEGER, are the order of the matrix A (N) and the
!          number of right-hand sides (M).
!
! D:       INTEGER, is the actual dimension of the arrays as they
!          are declared in the main program. D >= N is necessary.
!
! A:       REAL*8, declared as DIMENSION A(1:D,1:N), is the coefficient
!          matrix. On return it will be over-written by U.
!
! B:       REAL*8, declared as DIMENSION B(1:D,1:M), contains the right
!          hand sides (columnwise). On return each column occupies the
!          corresponding solution.
!
! EPS:     REAL*8, is the tolerance: if a pivotal element is
!          (absolutely) smaller than EPS, matrix A is considered
!          to be singular.
!
! DET:     REAL*8, will contain the determinant (DET = 0 is set,
!          if A is singular). It is scaled such that
!          1.0D-10 <= ABS(DET) <= 1.0D10 .
!
! EX:      INTEGER, is the scale factor for DET. The determinant
!          of A can be obtained from  DET*10**EX .
!
! CTR:     INTEGER, is a control variable with the following
!          meaning on input:
!          CTR < 0: evaluation of just the determinant. Parameter
!                   B is not used at all.
!          CTR = 0: solution of the system of equations and evaluation
!                   of the determinant.
!          CTR > 0: solution of the system of equations only. EX = 0
!                   and DET = 0 (1) will be set, if A is singular
!                   (regular).
!          On return it contains the error code:
!          CTR = -1: parameter fault (see REMARK 1).
!          CTR =  0: program worked.
!          CTR =  1: matrix A is singular.
!
! S:       INTEGER, declared as DIMENSION S(1:N), is used for
!          auxiliary storage only and has no meaning otherwise.
!
! REMARKS: (1) If N < 1 or M < 1, CTR = -1 is set and control is
!              returned to the calling program.
!          (2) The original contents of matrices A and B are
!              destroyed on return.
!          (3) Choosing B = I (unit matrix) with M = N will
!              generate the inverse matrix of A in B.
!          (4) This program makes no use of any special forms or
!              properties of the coefficient matrix (Hessenberg or
!              tridiagonal form, band matrices, symmetric or sparse
!              matrices). Special algorithms should be used in
!              that case, if the order is sufficiently high.
!
! REF.:    J.H. Wilkinson,
!          The Algebraic Eigenvalue Problem,
!          Oxford University Press, pp 212-214, 1965.
!
! AUTHOR:        H. RIEGER, RHRZ
! INSTALLATION:  IBM 3081, VS FORTRAN V2 (OPT=3)
!
! ACCESS:
!                --- MVS ---              --- VM/CMS ---
! LOAD MODULE:   SYS3.FORTLIB(DCOPIV)     FORTLIB  TXTLIB P (DCOPIV)
! DESCRIPTION:   SYS3.INFOLIB(DCOPIV)     INFOLIB  MACLIB P (DCOPIV)
! SOURCE MODULE: SYS3.SYMLIB.FORTRAN(DCOPIV)
!
!**********************************************************************

subroutine DCOPIV(A,B,N,M,iD,EPS,DET,iEX,iCTR,S)

implicit real*8(a-h,o-z)
dimension A(iD,N), B(iD,M), S(N)

call DCOPIV_INTERNAL(S)

! This is to allow type punning without an explicit interface
contains

subroutine DCOPIV_INTERNAL(S)

  use iso_c_binding

  real*8, target :: S(*)
  integer, pointer :: iS(:)

  call c_f_pointer(c_loc(S(1)),iS,[N])
  if ((N < 1) .or. (M < 1)) GO TO 7
  iEX = 0
  DET = 1.0d0
  if (N == 1) GO TO 4
  do I=1,N-1
    C = abs(A(I,I))
    iRI = I
    iCI = I
    do J=I,N
      do K=I,N
        if (abs(A(J,K)) > C) then
          C = abs(A(J,K))
          iRI = J
          iCI = K
        end if
      end do
    end do
    if (iRI /= I) then
      DET = -DET
      do J=1,N
        C = A(I,J)
        A(I,J) = A(iRI,J)
        A(iRI,J) = C
      end do
      if (iCTR >= 0) then
        do J=1,M
          C = B(I,J)
          B(I,J) = B(iRI,J)
          B(iRI,J) = C
        end do
      end if
    end if
    if (iCI /= I) then
      DET = -DET
      do J=1,N
        C = A(J,I)
        A(J,I) = A(J,iCI)
        A(J,iCI) = C
      end do
    end if
    iS(I) = iCI
    SUM = A(I,I)
    if (abs(SUM) <= EPS) then
      write(6,*) ' case 1. i,sum,eps',i,sum,eps
      GO TO 1
    end if
    do J=I+1,N
      C = A(J,I)/SUM
      do K=I+1,N
        A(J,K) = A(J,K)-C*A(I,K)
      end do
      if (iCTR >= 0) then
        do K=1,M
          B(J,K) = B(J,K)-C*B(I,K)
        end do
      end if
    end do
  end do
4 SUM = A(N,N)
  if (abs(SUM) <= EPS) then
    write(6,*) ' case 2. n,sum,eps',n,sum,eps
    GO TO 1
  end if
  if (iCTR > 0) GO TO 6
  do K=1,N
    DET = DET*A(K,K)
2   if (abs(DET) > 1.0d10) then
      DET = DET*1.0D-20
      iEX = iEX+20
      GO TO 2
    end if
3   if (abs(DET) <= 1.0D-10) then
      DET = DET*1.0d20
      iEX = iEX-20
      GO TO 3
    end if
  end do
  if (iCTR < 0) GO TO 5
6 do K=1,M
    B(N,K) = B(N,K)/SUM
  end do
  if (N == 1) GO TO 5
  do I=N-1,1,-1
    C = A(I,I)
    do J=1,M
      SUM = B(I,J)
      do K=I+1,N
        SUM = SUM-A(I,K)*B(K,J)
      end do
      B(I,J) = SUM/C
    end do
  end do
  do K=N-1,1,-1
    I = iS(K)
    if (I /= K) then
      do J=1,M
        C = B(I,J)
        B(I,J) = B(K,J)
        B(K,J) = C
      end do
    end if
  end do
5 iCTR = 0
  nullify(iS)
  return
1 DET = 0.0d0
  iCTR = 1
  nullify(iS)
  return
7 iCTR = -1
  nullify(iS)
  return
end subroutine DCOPIV_INTERNAL

end subroutine DCOPIV
