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
      SUBROUTINE HZLP2(CBUF,SBUF,DBUF,CSECT,RSECT,XI1,XI2,CNEW,ICI)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (ONE=1.0D00)

#include "SysDef.fh"

#include "mrci.fh"
      DIMENSION CBUF(MBUF,MXVEC),SBUF(MBUF,MXVEC),DBUF(MBUF),ICI(MBUF)
      DIMENSION CSECT(NSECT,MXVEC),RSECT(NSECT,MXVEC)
      DIMENSION CNEW(NSECT,MXVEC)
      DIMENSION XI1(NSECT,MXVEC),XI2(NSECT,MXVEC)
      DIMENSION IDCR(MXVEC),IDCW(MXVEC),IDS(MXVEC)
C THIS SUBROUTINE LOOPS OVER SECTIONS OF PSI AND SIGMA ARRAYS
C ON DISK, AND FORMS A NEW SET OF PSI ARRAYS AS A LINEAR COMBINATION
C OF THE BASIS SET PSI, RHO, XI1 AND XI2. TO FORM THE NEW PSI
C ARRAY, THE FIRST (NRROOT-NNEW) COLUMNS ARE SKIPPED, SINCE THEY
C GIVE NO ESSENTIAL IMPROVEMENT.
      IVZSTA=1+NRROOT-NNEW
      IVZ1=1
      IVZ2=1+NRROOT
      IVZ3=1+2*NRROOT
      IVZ4=1+3*NRROOT
C     WRITE(6,*)
C     WRITE(6,*)' IN HZLP2. NNEW=',NNEW
C     IF(NVEC.LT.MXVEC) WRITE(6,*)' DUMMY WRITES OF NEW FUNCTIONS.'
C WE MAY NEED DUMMY WRITES TO PROVIDE DISK ADDRESSES:
      NNVEC=MIN(NVEC+NNEW,MXVEC)
      DO 16 K=NVEC+1,NNVEC
        IDISKC(K)=IDFREE
C       WRITE(6,'(A,I2,A,I8)')' IDISKC(',K,')=',IDFREE
        DO 15 ISTA=1,NCONF,MBUF
          IEND=MIN(NCONF,ISTA+MBUF-1)
          IBUF=1+IEND-ISTA
          CALL iDAFILE(LUEIG,0,ICI,IBUF,IDFREE)
15      CONTINUE
16    CONTINUE
C WE NEED COPIES OF THE DISK ADDRESSES. TWO COPIES FOR PSI BUFFERS.
      DO 10 K=1,NNVEC
        IDCR(K)= IDISKC(K)
        IDCW(K)= IDISKC(K)
        IDS(K)= IDISKS(K)
10    CONTINUE
      IDD=IDISKD
C LOOP OVER BUFFERS FOR READING PSI, SIGMA AND DBUF:
      DO 2000 ISTA=1,NCONF,MBUF
        IEND=MIN(NCONF,ISTA+MBUF-1)
        IBUF=1+IEND-ISTA
        CALL dDAFILE(LUEIG,2,DBUF,IBUF,IDD)
        DO 20 K=1,NVEC
          CALL iDAFILE(LUEIG,2,ICI,IBUF,IDCR(K))
          CALL UPKVEC(IBUF,ICI,CBUF(1,K))
          CALL dDAFILE(LUEIG,2,SBUF(1,K),IBUF,IDS(K))
20      CONTINUE
C LOOP OVER VECTOR SECTIONS, LENGTH AT MOST NSECT:
        DO 1000 JSTA=1,IBUF,NSECT
          JEND=MIN(IBUF,JSTA+NSECT-1)
          ISECT=1+JEND-JSTA
C TRANSFORM TO EIGENFUNCTIONS OF HSMALL: FIRST, CI SECTION.
        CALL DGEMM_('N','N',
     &              ISECT,NRROOT,NVEC,
     &              1.0d0,CBUF(JSTA,1),MBUF,
     &              VSMALL,MXVEC,
     &              0.0d0,CSECT,NSECT)
C THEN, SIGMA SECTION INTO RSECT.
        CALL DGEMM_('N','N',
     &              ISECT,NRROOT,NVEC,
     &              1.0d0,SBUF(JSTA,1),MBUF,
     &              VSMALL,MXVEC,
     &              0.0d0,RSECT,NSECT)
C AND THEN FORM RSECT=SECTION OF RESIDUAL ARRAY, AND XI1 AND XI2:
        DO 30 I=1,ISECT
          DO 31 K=1,NRROOT
            RSECT(I,K)=RSECT(I,K)-ESMALL(K)*CSECT(I,K)
            XI1(I,K)=CSECT(I,K)/(DBUF(I+JSTA-1)-ESMALL(K))
            XI2(I,K)=RSECT(I,K)/(DBUF(I+JSTA-1)-ESMALL(K))
31        CONTINUE
30      CONTINUE
C FORM NEW PSI ARRAYS IN CNEW SECTION:
          CALL DGEMM_('N','N',
     &                ISECT,NNEW,NRROOT,
     &                1.0d0,CSECT,NSECT,
     &                VZERO(IVZ1,IVZSTA),MXZ,
     &                0.0d0,CNEW,NSECT)
          CALL DGEMM_('N','N',ISECT,NNEW,NRROOT,ONE,RSECT,NSECT,
     *          VZERO(IVZ2,IVZSTA),MXZ,ONE,CNEW,NSECT)
          CALL DGEMM_('N','N',ISECT,NNEW,NRROOT,ONE,XI1  ,NSECT,
     *          VZERO(IVZ3,IVZSTA),MXZ,ONE,CNEW,NSECT)
          CALL DGEMM_('N','N',ISECT,NNEW,NRROOT,ONE,XI2  ,NSECT,
     *          VZERO(IVZ4,IVZSTA),MXZ,ONE,CNEW,NSECT)
C     IF(ISTA+JSTA.EQ.2) THEN
C       WRITE(6,*)' CONSTRUCTION OF NEW VECTOR IN HZLP2.'
C       WRITE(6,*)' CSECT:'
C       WRITE(6,'(1X,5F15.6)')((CSECT(I,J),I=1,5),J=1,NNEW)
C       WRITE(6,*)' RSECT:'
C       WRITE(6,'(1X,5F15.6)')((RSECT(I,J),I=1,5),J=1,NNEW)
C       WRITE(6,*)'   XI1:'
C       WRITE(6,'(1X,5F15.6)')((XI1(I,J),I=1,5),J=1,NNEW)
C       WRITE(6,*)'   XI2:'
C       WRITE(6,'(1X,5F15.6)')((XI2(I,J),I=1,5),J=1,NNEW)
C       WRITE(6,*)' VZERO:'
C       IIII=IVZSTA-1
C       WRITE(6,'(1X,4F15.6)')((VZERO(I,IIII+J),I=1,4*NNEW),J=1,NNEW)
C       WRITE(6,*)'  CNEW:'
C       WRITE(6,'(1X,5F15.6)')((CNEW(I,J),I=1,5),J=1,NNEW)
C     END IF
C INSERT THE NEW PSI SECTIONS IN BUFFER. THIS MAY IMPLY OVERWRITING
C OLD ENTRIES, BUT CAN ALSO LEAD TO AN INCREASED NUMBER OF VECTORS:
          DO 50 KK=1,NNEW
            NN=NVTOT+KK
            K=1+MOD(NN-1,MXVEC)
C         IF(ISTA+JSTA.EQ.2) THEN
C           WRITE(6,'(A,I2,A,I6)')' CNEW NR.',KK,' COPIED TO BUFFER ',K
C           WRITE(6,*)' IT CONTAINS:'
C           WRITE(6,'(1X,5F15.6)')(CNEW(I,KK),I=1,15)
C         END IF
            CALL DCOPY_(ISECT,CNEW(1,KK),1,CBUF(JSTA,K),1)
50        CONTINUE
C CONTINUE, NEXT SECTION.
1000    CONTINUE
        DO 60 KK=1,NNEW
          NN=NVTOT+KK
          K=1+MOD(NN-1,MXVEC)
C         IF(ISTA.EQ.1) THEN
C           WRITE(6,'(A,I2,A,I6)')' BUFFER NR.',K,' WRITTEN AT ',IDCW(K)
C           WRITE(6,*)' IT CONTAINS:'
C           WRITE(6,'(1X,5F15.6)')(CBUF(I,K),I=1,15)
C         END IF
          CALL PKVEC(IBUF,CBUF(1,K),ICI)
          CALL iDAFILE(LUEIG,1,ICI,IBUF,IDCW(K))
60      CONTINUE
C CONTINUE, NEXT BUFFER.
2000  CONTINUE
      NVTOT=NVTOT+NNEW
      RETURN
      END
