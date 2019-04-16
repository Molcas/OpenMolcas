************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      SUBROUTINE SECULAR(NDIM,N,NRON,HMAT,SMAT,VEC,EVAL,SCR,THR)
      IMPLICIT REAL*8 (A-H,O-Z)
      INTRINSIC SQRT
      DIMENSION HMAT(NDIM,NDIM),SMAT(NDIM,NDIM)
      DIMENSION VEC(NDIM,NDIM),EVAL(NDIM),SCR(*)
      THR2=THR**2
C PUT NORMALIZED VECTORS INTO VEC:
      CALL DCOPY_(N*NDIM,[0.0D00],0,VEC,1)
      DO 5 I=1,N
        VEC(I,I)=1.0D00/SQRT(SMAT(I,I))
5     CONTINUE
C GRAM-SCHMIDT ORTHONORMALIZING PROCEDURE:
      NRON=0
      DO 60 I=1,N
C SMAT*(NORMALIZED VECTOR) INTO SCR:
        CALL DCOPY_(N,SMAT(1,I),1,SCR,1)
        CALL DSCAL_(N,VEC(I,I),SCR,1)
C PROJECT AWAY THE ALREADY ORTHONORMALIZED BASIS SET:
        DO 30 J=1,NRON
          MAXLEN=I-1-NRON+J
          SUM=0.0D00
          DO 10 K=1,MAXLEN
            SUM=SUM+VEC(K,J)*SCR(K)
10        CONTINUE
          DO 20 K=1,MAXLEN
            VEC(K,I)=VEC(K,I)-SUM*VEC(K,J)
20        CONTINUE
30      CONTINUE
C NORMALIZE AND MOVE INTO POSITION:
        SUM=0.0D00
        DO 40 K=1,I
          SUM=SUM+VEC(K,I)*SCR(K)
40      CONTINUE
        IF(SUM.LT.THR2) GOTO 60
        NRON=NRON+1
        SCALE=1.0D00/SQRT(SUM)
        DO 50 K=1,I
          VEC(K,NRON)=SCALE*VEC(K,I)
50      CONTINUE
60    CONTINUE
      DO 70 I=NRON+1,N
        CALL DCOPY_(N,[0.0D00],0,VEC(1,I),1)
70    CONTINUE
C TRANSFORM HAMILTONIAN INTO SCR:
      IOFF1=N*NRON
      CALL DGEMM_('N','N',
     &            N,NRON,N,
     &            1.0d0,HMAT,NDIM,
     &            VEC,NDIM,
     &            0.0d0,SCR,N)
      CALL DGEMM_('T','N',
     &            NRON,NRON,N,
     &            1.0d0,VEC,NDIM,
     &            SCR,N,
     &            0.0d0,SCR(IOFF1+1),NRON)
C COPY TRANSFORMED HMAT INTO TRIANGULAR STORAGE IN SCR:
      IFROM=IOFF1+1
      ITO=1
      DO 110 I=1,NRON
        CALL DCOPY_(I,SCR(IFROM),1,SCR(ITO),1)
        ITO=ITO+I
        IFROM=IFROM+NRON
110   CONTINUE
C DIAGONALIZE:
      CALL Jacob(SCR,VEC,NRON,NDIM)
C COPY EIGENVALUES INTO EVAL:
      II=0
      DO 120 I=1,NRON
        II=II+I
        EVAL(I)=SCR(II)
120   CONTINUE
      RETURN
      END
