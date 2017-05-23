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
      SUBROUTINE HZLP1(CBUF,SBUF,DBUF,ARR,CSECT,RSECT,XI1,XI2,ICI)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (ONE=1.0D00)
      PARAMETER (IX1F=1,IX2F=2,IRR=3,IX1R=4,IX2R=5,IX1X1=6,IX2X1=7,
     *           IX2X2=8,IFDF=9,IFDR=10,IRDR=11)

#include "SysDef.fh"

#include "mrci.fh"
      DIMENSION CBUF(MBUF,MXVEC),SBUF(MBUF,MXVEC),DBUF(MBUF),ICI(MBUF)
      DIMENSION CSECT(NSECT,MXVEC),RSECT(NSECT,MXVEC)
      DIMENSION XI1(NSECT,NRROOT),XI2(NSECT,NRROOT)
      DIMENSION ARR(NRROOT,NRROOT,11)
      DIMENSION IDCOPY(MXVEC),IDC(MXVEC),IDS(MXVEC)
C THIS SUBROUTINE LOOPS OVER SECTIONS OF PSI AND SIGMA ARRAYS
C ON DISK, AND ACCUMULATES OVERLAP MATRICES AND A COUPLE OF
C HAMILTONIAN MATRICES IN THE BASIS SET PSI, RHO, XI1 AND XI2.
C THE 11 MATRICES X1F,..,RDR ARE STORED CONSECUTIVELY IN THE
C SINGLE ARRAY ARR.
      NRR2=NRROOT**2
      CALL DCOPY_(11*NRR2,0.0D00,0,ARR,1)
      DO 10 K=1,NVEC
        IDC(K)= IDISKC(K)
        IDS(K)= IDISKS(K)
10    CONTINUE
      IDD= IDISKD
C LOOP OVER BUFFERS FOR READING PSI, SIGMA AND DBUF:
      DO 2000 ISTA=1,NCONF,MBUF
        IEND=MIN(NCONF,ISTA+MBUF-1)
        IBUF=1+IEND-ISTA
        CALL dDAFILE(LUEIG,2,DBUF,IBUF,IDD)
        DO 20 K=1,NVEC
          IDCOPY(K)=IDC(K)
          CALL iDAFILE(LUEIG,2,ICI,IBUF,IDC(K))
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
          DO 30 K=1,NRROOT
            RSECT(I,K)=RSECT(I,K)-ESMALL(K)*CSECT(I,K)
            XI1(I,K)=CSECT(I,K)/(DBUF(I+JSTA-1)-ESMALL(K))
            XI2(I,K)=RSECT(I,K)/(DBUF(I+JSTA-1)-ESMALL(K))
30      CONTINUE
C ACCUMULATE OVERLAP MATRICES:
        CALL DGEMM_('T','N',NRROOT,NRROOT,ISECT,ONE,XI1,NSECT,
     *           CSECT,NSECT,ONE,ARR(1,1,IX1F),NRROOT)
        CALL DGEMM_('T','N',NRROOT,NRROOT,ISECT,ONE,XI2,NSECT,
     *           CSECT,NSECT,ONE,ARR(1,1,IX2F),NRROOT)
        CALL DGEMM_('T','N',NRROOT,NRROOT,ISECT,ONE,RSECT,NSECT,
     *           RSECT,NSECT,ONE,ARR(1,1,IRR),NRROOT)
        CALL DGEMM_('T','N',NRROOT,NRROOT,ISECT,ONE,XI1,NSECT,
     *           RSECT,NSECT,ONE,ARR(1,1,IX1R),NRROOT)
        CALL DGEMM_('T','N',NRROOT,NRROOT,ISECT,ONE,XI2,NSECT,
     *           RSECT,NSECT,ONE,ARR(1,1,IX2R),NRROOT)
        CALL DGEMM_('T','N',NRROOT,NRROOT,ISECT,ONE,XI1,NSECT,
     *           XI1  ,NSECT,ONE,ARR(1,1,IX1X1),NRROOT)
        CALL DGEMM_('T','N',NRROOT,NRROOT,ISECT,ONE,XI2,NSECT,
     *           XI1  ,NSECT,ONE,ARR(1,1,IX2X1),NRROOT)
        CALL DGEMM_('T','N',NRROOT,NRROOT,ISECT,ONE,XI2,NSECT,
     *           XI2  ,NSECT,ONE,ARR(1,1,IX2X2),NRROOT)
C PUT D*CSECT INTO XI1, AND D*RSECT INTO XI2:
        DO 40 I=1,ISECT
          DO 40 K=1,NRROOT
            XI1(I,K)=DBUF(I+JSTA-1)*CSECT(I,K)
            XI2(I,K)=DBUF(I+JSTA-1)*RSECT(I,K)
40      CONTINUE
C ACCUMULATE ARRAYS FDF, FDR, AND RDR:
        CALL DGEMM_('T','N',NRROOT,NRROOT,ISECT,ONE,XI1,NSECT,
     *           CSECT,NSECT,ONE,ARR(1,1,IFDF),NRROOT)
        CALL DGEMM_('T','N',NRROOT,NRROOT,ISECT,ONE,XI1,NSECT,
     *           RSECT,NSECT,ONE,ARR(1,1,IFDR),NRROOT)
        CALL DGEMM_('T','N',NRROOT,NRROOT,ISECT,ONE,XI2,NSECT,
     *           RSECT,NSECT,ONE,ARR(1,1,IRDR),NRROOT)
C CONTINUE, NEXT SECTION.
1000    CONTINUE
C CONTINUE, NEXT BUFFER.
2000  CONTINUE
      RETURN
      END
