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
      SUBROUTINE CPLAB (A,B,L,M,N,IA,IB,C,IC,IER)
      implicit real*8(a-h,o-z)
!ulf      INTEGER            L,M,N,IA,IB,IC,IER
!ulf      real*8   A(IA,M),B(IB,N),C(IC,N)
!ulf      real*8   TEMP
      DIMENSION   A(IA,M),B(IB,N),C(IC,N)
      IF (IA .GE. L .AND. IB .GE. M .AND. IC .GE. L) GO TO 5
      IER=129
      GO TO 9000
    5 IER = 0
      DO 15 I=1,L
         DO 16 J=1,N
            TEMP=0.0D0
            DO 10 K=1,M
               TEMP=A(I,K)*B(K,J)+TEMP
   10       CONTINUE
            C(I,J)=C(I,J)+TEMP
   16    CONTINUE
   15 CONTINUE
      GO TO 9005
 9000 CONTINUE
 9005 RETURN
      END
