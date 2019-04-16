************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 1986, Per E. M. Siegbahn                               *
*               1986, Margareta R. A. Blomberg                         *
************************************************************************
      SUBROUTINE FMUL2_CPF(A,B,C,NROW,NCOL,N)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(NROW,N),B(NCOL,N),CJ(200)
      DIMENSION C(NROW,NCOL)

      If ( nRow.gt.200 ) then
         WRITE(6,*)
      CALL XFLUSH(6)
         WRITE(6,*) ' *** Error in Subroutine FMUL2_CPF ***'
      CALL XFLUSH(6)
         WRITE(6,*) ' row dimension exceeds local buffer size'
      CALL XFLUSH(6)
         WRITE(6,*)
      CALL XFLUSH(6)
         CALL Abend
      End If

      DO 10 J=1,NCOL
      DO 15 I=1,NROW
      CJ(I)=0.0D0
15    CONTINUE
      IF(J.EQ.NCOL)GO TO 16
      J1=J+1
      DO 20 K=1,N
      FAC=B(J,K)
      IF(FAC.EQ.0.0)GO TO 20
      DO 25 I=J1,NROW
      CJ(I)=CJ(I)+FAC*A(I,K)
25    CONTINUE
20    CONTINUE
16    DO 30 I=1,NROW
      C(I,J)=CJ(I)
30    CONTINUE
10    CONTINUE
      RETURN
      END
