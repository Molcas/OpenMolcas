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
      SUBROUTINE MYDGEMM ( DoIt, M, N, K,                               &
     &                     A, LDA, B, LDB,                              &
     &                     C, LDC )
      Use Constants, only: Zero
      Implicit None
!     .. Scalar Arguments ..
      INTEGER            M, N, K, LDA, LDB, LDC
      Integer DoIt(*)
!     .. Array Arguments ..
      REAL*8   A( LDA, * ), B( LDB, * ), C( LDC, * )
!     ..
!  Purpose
!  =======
!
!  DGEMM for a special case.
!
!
!     .. Local Scalars ..
      INTEGER            J, L
!
!     Form  C := A*B + C.
!
      DO J = 1, N
         If (DoIt(J).ne.1) Cycle
         DO L = 1, K
            IF( B( L, J ).EQ.ZERO ) Cycle
            Call DAxPy_(M,B(L,J),A(:,L),1,C(:,J),1)
         End Do
      End Do
!
      END SUBROUTINE MYDGEMM
