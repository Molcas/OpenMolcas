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
      SUBROUTINE ZJAC(NDIMEN,ARRRE,ARRIM,LDV,VECRE,VECIM)
      IMPLICIT NONE
      INTEGER NDIMEN,LDV
      REAL*8 ARRRE(NDIMEN,NDIMEN),ARRIM(NDIMEN,NDIMEN)
      REAL*8 VECRE(LDV,*),VECIM(LDV,*)
      REAL*8 EPS
      PARAMETER (EPS=1.0D-12)
      REAL*8 SBDMAX,ARII,ARJJ,ARIJ,AIIJ,AAIJ,ERE,EIM
      REAL*8 DIFF,SGN,DUM,CS,SN,TN
      REAL*8 VRKJ,VIKJ,VRKI,VIKI,ARKJ,AIKJ,ARKI,AIKI
      REAL*8 ARJK,AIJK,ARIK,AIIK,ESEL,EVAL
      REAL*8 VNSUM,VNOLD
      INTEGER I,J,K,ISEL
      INTEGER NR,NROT,NSWEEP,IFTEST,IFERR
      REAL*8 ERRRE,ERRIM
      REAL*8 O_i, O_j, Thr_EDiff, EDIFF

      IFTEST=0

C von Neumanns sum should be ever decreasing. Check this:
      VNSUM=1.0D99
C Max sub-diagonal element:
      NROT=0
      NSWEEP=0
  10  CONTINUE
      NSWEEP=NSWEEP+1
      NR=0
      SBDMAX=EPS
      VNOLD=VNSUM
      VNSUM=0.0D0
      DO I=2,NDIMEN
       DO J=1,I-1
        ARII=ARRRE(I,I)
        ARJJ=ARRRE(J,J)
        ARIJ=ARRRE(I,J)
        AIIJ=ARRIM(I,J)
        AAIJ=SQRT(ARIJ**2+AIIJ**2)
        VNSUM=VNSUM+AAIJ**2
        SBDMAX=MAX(SBDMAX,AAIJ)
        IF(2.0D0*AAIJ.LT.SBDMAX) GOTO 100
        NR=NR+1
        ERE=ARIJ/AAIJ
        EIM=AIIJ/AAIJ
        DIFF=ARII-ARJJ
        SGN=1.0D0
        IF(DIFF.LT.0.0D0) THEN
           DIFF=-DIFF
           SGN=-SGN
        END IF
        DUM=DIFF+SQRT(DIFF**2+4.0D0*AAIJ**2)
        TN=2.0D0*SGN*AAIJ/DUM
        CS=1.0D0/SQRT(1.0D0+TN**2)
        SN=CS*TN
C TN,CS,SN=TAN,COS AND SIN OF ROTATION ANGLE.
C Rotate vectors:
        DO K=1,LDV
         VRKJ=CS*VECRE(K,J)-SN*( ERE*VECRE(K,I)-EIM*VECIM(K,I))
         VIKJ=CS*VECIM(K,J)-SN*( EIM*VECRE(K,I)+ERE*VECIM(K,I))
         VRKI=SN*VECRE(K,J)+CS*( ERE*VECRE(K,I)-EIM*VECIM(K,I))
         VIKI=SN*VECIM(K,J)+CS*( EIM*VECRE(K,I)+ERE*VECIM(K,I))
         VECRE(K,J)=VRKJ
         VECIM(K,J)=VIKJ
         VECRE(K,I)=VRKI
         VECIM(K,I)=VIKI
        END DO
C Rotate matrix columns:
        DO K=1,NDIMEN
         ARKJ=CS*ARRRE(K,J)-SN*( ERE*ARRRE(K,I)-EIM*ARRIM(K,I))
         AIKJ=CS*ARRIM(K,J)-SN*( EIM*ARRRE(K,I)+ERE*ARRIM(K,I))
         ARKI=SN*ARRRE(K,J)+CS*( ERE*ARRRE(K,I)-EIM*ARRIM(K,I))
         AIKI=SN*ARRIM(K,J)+CS*( EIM*ARRRE(K,I)+ERE*ARRIM(K,I))
         ARRRE(K,J)=ARKJ
         ARRIM(K,J)=AIKJ
         ARRRE(K,I)=ARKI
         ARRIM(K,I)=AIKI
        END DO
C Rotate matrix rows:
        DO K=1,NDIMEN
         ARJK=CS*ARRRE(J,K)-SN*( ERE*ARRRE(I,K)+EIM*ARRIM(I,K))
         AIJK=CS*ARRIM(J,K)-SN*(-EIM*ARRRE(I,K)+ERE*ARRIM(I,K))
         ARIK=SN*ARRRE(J,K)+CS*( ERE*ARRRE(I,K)+EIM*ARRIM(I,K))
         AIIK=SN*ARRIM(J,K)+CS*(-EIM*ARRRE(I,K)+ERE*ARRIM(I,K))
         ARRRE(J,K)=ARJK
         ARRIM(J,K)=AIJK
         ARRRE(I,K)=ARIK
         ARRIM(I,K)=AIIK
        END DO
 100    CONTINUE
       END DO
      END DO
      NROT=NROT+NR

      IF(IFTEST.GT.0) THEN
C --- CHECK IF DIVERGING (This should never happen):
       IF(VNSUM.GE.VNOLD) THEN
        Call WarningMessage(2,'ZJAC got increasing von-Neumann-sum.')
        WRITE(6,*)' Panic exit from Jacobi iteration loop.'
        GOTO 999
       END IF
C --- CHECK IF IDLE LOOPS (This should never happen):
       IF(NR.EQ.0 .AND. SBDMAX.GT.EPS) THEN
        Call WarningMessage(2,'ZJAC detected infinite idle loops.')
        WRITE(6,*)' Panic exit from Jacobi iteration loop.'
        GOTO 999
       END IF
      END IF

* PAM 2008 If unchecked, then idle loop just generates warning and return:
      IFERR=0
      IF(IFTEST.EQ.0) THEN
       IF(VNSUM.GE.VNOLD) THEN
        Call WarningMessage(1,'ZJAC got increasing von-Neumann-sum.')
        IFERR=1
       END IF
       IF(NR.EQ.0 .AND. SBDMAX.GT.EPS) THEN
        Call WarningMessage(1,'ZJAC detected infinite idle loops.')
        IFERR=1
       END IF
       IF(IFERR.NE.0) THEN
        WRITE(6,*)' Probably, the convergence criteria of the ZJAC'
        WRITE(6,*)' code are too strict. You may report this as a bug'
        WRITE(6,*)' but since you probably simply hit the limits to'
        WRITE(6,*)' accuracy caused by round-off, the program will'
        WRITE(6,*)' continue.'
        WRITE(6,*)' At this point, the largest subdiagonal element'
        WRITE(6,*)' is SBDMAX=',SBDMAX
       END IF
      END IF

C --- CHECK IF CONVERGED:
      IF(IFERR.EQ.0 .AND. SBDMAX.GT.EPS) GOTO 10

C Order the eigenvalues in increasing sequence:
C
C     Note that special care is taken to define some order in case
C     of degeneracy.
C
      Thr_EDiff=1.0D-10
      DO I=1,NDIMEN-1
         ESEL=ARRRE(I,I)
         ISEL=I
*
         O_i=0.0D0
         Do K=1,LDV
            O_i = O_i + DBLE(K) * (VECRE(K,I)**2+VECIM(K,I)**2)
         End Do
*
*        Check against other eigenvalues
*
         DO J=I+1,NDIMEN
            EVAL=ARRRE(J,J)
*
            EDIFF=ABS(EVAL-ESEL)
            IF (EVAL.LT.ESEL.and.EDIFF.gt.Thr_EDiff) THEN
               ISEL=J
               ESEL=EVAL
            Else If (EDIFF.lt.Thr_EDiff) THEN
               O_j=0.0D0
               Do K=1,LDV
                  O_j = O_j + DBLE(K) * (VECRE(K,J)**2+VECIM(K,J)**2)
               End Do
               If (O_j.gt.O_i) Then
                  ISEL=J
                  ESEL=EVAL
               End If
            END IF
         END DO
*
*        Swap here!
*
         IF (ISEL.NE.I) THEN
            DO K=1,LDV
               VRKI=VECRE(K,I)
               VRKJ=VECRE(K,ISEL)
               VIKI=VECIM(K,I)
               VIKJ=VECIM(K,ISEL)
               VECRE(K,I)=VRKJ
               VECRE(K,ISEL)=VRKI
               VECIM(K,I)=VIKJ
               VECIM(K,ISEL)=VIKI
            END DO
            ARRRE(ISEL,ISEL)=ARRRE(I,I)
            ARRRE(I,I)=ESEL
         END IF
      END DO
*
      RETURN
 999  CONTINUE
C Jump here on error.
      Call WarningMessage(2,'ZJAC abend diagnostics:')
      WRITE(6,*)' Nr of sweeps      NSWEEP=',NSWEEP
      WRITE(6,*)' Nr of 2-rotations   NROT=',NROT
      WRITE(6,*)' Rotations this sweep  NR=',NR
      WRITE(6,*)' von Neumann sum    VNSUM=',VNSUM
      WRITE(6,*)' Previous value     VNOLD=',VNOLD
      WRITE(6,*)' Largest subdiag   SBDMAX=',SBDMAX
      CALL HRMCHK(NDIMEN,ARRRE,ARRIM,ERRRE,ERRIM)
      WRITE(6,*)' Antisymm. part of ARRRE:',ERRRE
      WRITE(6,*)' Symmetric part of ARRIM:',ERRIM
      CALL ABEND()
      END
      SUBROUTINE HRMCHK(NDIMEN,ARRRE,ARRIM,ERRRE,ERRIM)
      IMPLICIT NONE
      INTEGER NDIMEN,I,J
      REAL*8 ARRRE(NDIMEN,NDIMEN),ARRIM(NDIMEN,NDIMEN)
      REAL*8 ERRRE,ERRIM
      ERRRE=0.0D0
      ERRIM=0.0D0
      DO I=1,NDIMEN
       DO J=1,I
        ERRRE=MAX(ERRRE,ABS(ARRRE(I,J)-ARRRE(J,I)))
        ERRIM=MAX(ERRIM,ABS(ARRIM(I,J)+ARRIM(J,I)))
       END DO
      END DO
      RETURN
      END
