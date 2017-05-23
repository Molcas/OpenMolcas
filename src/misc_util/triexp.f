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
* Copyright (C) 1989, Markus P. Fuelscher                              *
************************************************************************
      SUBROUTINE TRIEXP(A,B,NDIM)
C
C     AUTHOR:        M.P. FUELSCHER, UNIV. OF LUND, SWEDEN, 1989
C     MODIFICATIONS: INCLUSION INTO THE RASSCF METHOD
C                    M.P. FUELSCHER, UNIV. OF LUND, SWEDEN, MAY 1989
C
C     EXPAND A REAL, SYMMETRIC MATRIX GIVEN AS THE LOWER TRIANGLE
C     INTO FULL STORAGE MODE.
C
C     CALLING PARAMETERS:
C     A   : INPUT MATRIX (LOWER TRINAGULAR STORAGE MODE)
C     B   : OUTPUT MATRIX (FULL STORAGE MODE)
C     NDIM: DIMENSION OF MATRIX A AND B
C
      REAL*8 A(*)
      REAL*8 B(*)
C
      IMODE=1
      If ( idLOC(A(1)).EQ.idLOC(B(1)) ) IMODE=2
C
      IF( IMODE.EQ.1 ) THEN
         K=0
         DO 10 I=1,NDIM
            DO 20 J=1,I
               K=K+1
               L1=J+(I-1)*NDIM
               B(L1)=A(K)
               L2=I+(J-1)*NDIM
               B(L2)=A(K)
20          CONTINUE
10       CONTINUE
      ENDIF
C
      IF( IMODE.EQ.2 ) THEN
         K=1+NDIM*(NDIM+1)/2
         DO 30 I=NDIM,1,-1
            DO 40 J=I,1,-1
               K=K-1
               L1=J+(I-1)*NDIM
               L2=I+(J-1)*NDIM
               L3=MAX(L2,L1)
               B(L3)=A(K)
40          CONTINUE
30       CONTINUE
         DO 50 I=1,NDIM
            DO 60 J=1,I
               L1=J+(I-1)*NDIM
               L2=I+(J-1)*NDIM
               L3=MAX(L2,L1)
               L4=MIN(L2,L1)
               B(L4)=B(L3)
60          CONTINUE
50       CONTINUE
      ENDIF
C
      RETURN
      END
