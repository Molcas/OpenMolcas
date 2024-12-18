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
* Copyright (C) 2000, Per Ake Malmqvist                                *
************************************************************************
*****************************************************************
*  PROGRAM RASSI        PER-AAKE MALMQVIST 2000-06-30
*  SUBROUTINE TRAINT
*  TRANSFORM TWO-ELECTRON INTEGRALS TO A MIXED ACTIVE ORBITAL
*  BASIS. THE RESULTING INTEGRALS ARE STORED IN ARRAY TUVX IN
*  THE SAME FORMAT AS THE 2-EL TRANSITION DENSITY MATRICES.
*****************************************************************
      SUBROUTINE TRAINT(CMO1,CMO2,NGAM2,TUVX)
      use Constants, only: Zero
      use stdalloc, only: mma_allocate, mma_deallocate
      use TRNSFRM, only: NAP,NBP,ISP,LMOP1,NAQ,NBQ,ISQ,LMOQ1,NAR,NBR,
     &                   LMOR1,ISR,NAS,NBS,LMOS1,ISS,NBPQ,NBRS,NAVX,
     &                   NX1MX,NX2MX,NX3MX,NVXPQ,IAPR
      use Symmetry_Info, only: nSym=>nIrrep, MUL
      use rassi_data, only: NCMO,NBMX,NASH,NBASF,NISH,NOSH
      IMPLICIT NONE
      INTEGER NGAM2
      Real*8 CMO1(NCMO),CMO2(NCMO),TUVX(NGAM2)

      INTEGER KEEP(8),NBSX(8)
      LOGICAL   ISQARX
      Real*8, Allocatable:: X1(:), X2(:), X3(:), VXPQ(:)
      INTEGER IRC,IA,IS,INTBUF,LMOP,NSP,LMOQ,NSQ,NSPQ,ISPQ,LMOR,NSRM,
     &        NSR,NSPQR,LMOS,NSSM,NSS,ISRS,NACT,NSYMX
C CLEAR THE ARRAY OF TRANSFORMED INTEGRALS.
      TUVX(:)=Zero
C RETRIEVE STRUCTURE DATA FOR THE ORDERED INTEGRAL FILE.
C KEEP(IS) IS POSITIVE TO INDICATE THAT BASIS FUNCTIONS WITH
C SYMMETRY LABEL IS SHOULD BE SKIPPED. KEEP(2)=0 OR 1 TO SHOW
C IF SYMMETRY LABELS IJKL OF EACH SYMMETRY BLOCK ARE CANONICAL
C (ONLY CASES I>=J,K>=L,IJ>=KL PRESENT) OR IF LEFT AND RIGHT HAND
C PAIR ARE INDIVIDUALLY TREATED (I>=J,K>=L), RESPECTIVELY.
C THE IDATA ARRAY GIVES DISK ADDRESS, AND NUMBER OF INTEGRAL MATRICES
C PER BUFFER, FOR EACH SYMMETRY BLOCK. THERE ARE AT MOST 176 SUCH
C BLOCKS, AND THEIR ORDERING IS DETERMINED BY THE SYMMETRY
C LOOPS, WHICH MUST BE THE SAME AS IN THE ORDERING PROGRAM.
C RETRIEVE BASE DATA FROM UNIT LUORD:
C RETRIEVE BASE DATA FROM UNIT LUORD:
      IRC=0
      CALL GETORD(IRC,ISQARX,NSYMX,NBSX,KEEP)
C SET UP IAPR(IS)=NR OF ACTIVES WITH PREVIOUS SYMMETRY LABEL.
      IA=0
      DO 10 IS=1,NSYM
        IAPR(IS)=IA
        IA=IA+NASH(IS)
10    CONTINUE
C LOOP OVER QUADRUPLES OF SYMMETRIES (NSP,NSQ,NSR,NSS) NSR>=NSS
C IN THE SAME ORDER AS IN THE INTORD PROGRAM.
      INTBUF=MAX(NBMX**2,256*256)
      LMOP=1
      DO 104 NSP=1,NSYM
       IF(NSP.NE.1) LMOP=LMOP+NBASF(NSP-1)*NOSH(NSP-1)
       NAP=NASH(NSP)
C      KEEPP=KEEP(NSP)
       NBP=NBASF(NSP)
       ISP=NSP
       LMOP1=LMOP+NISH(NSP)*NBP
       LMOQ=1
       DO 103 NSQ=1,NSP
        IF(NSQ.NE.1) LMOQ=LMOQ+NBASF(NSQ-1)*NOSH(NSQ-1)
        NAQ=NASH(NSQ)
C       KEEPQ=KEEP(NSQ)
        NBQ=NBASF(NSQ)
        NSPQ=MUL(NSP,NSQ)
        ISQ=NSQ
        ISPQ=(ISP**2-ISP)/2+ISQ
        LMOQ1=LMOQ+NISH(NSQ)*NBQ
        LMOR=1
        NSRM=NSYM
        IF(ISQARX) NSRM=NSP
        DO 102 NSR=1,NSRM
         IF(NSR.NE.1) LMOR=LMOR+NBASF(NSR-1)*NOSH(NSR-1)
         NAR=NASH(NSR)
C        KEEPR=KEEP(NSR)
         NBR=NBASF(NSR)
         LMOR1=LMOR+NBR*NISH(NSR)
         NSPQR=MUL(NSPQ,NSR)
         ISR=NSR
         LMOS=1
         NSSM=NSR
         DO 101 NSS=1,NSSM
          IF(NSS.NE.1) LMOS=LMOS+NBASF(NSS-1)*NOSH(NSS-1)
          IF(NSPQR.NE.NSS) GO TO 101
          NAS=NASH(NSS)
C         KEEPS=KEEP(NSS)
          NBS=NBASF(NSS)
          LMOS1=LMOS+NBS*NISH(NSS)
          ISS=NSS
          ISRS=(ISR**2-ISR)/2+ISS
C SHOULD THIS SYMMETRY BLOCK BE USED...?
C         KEEPT=KEEPP+KEEPQ+KEEPR+KEEPS
          NACT=NAP*NAQ*NAR*NAS
C ...WELL, IT IS PRESENT ON THE FILE
          IF(ISPQ.LT.ISRS) GO TO 101
          IF(NACT.EQ.0) GO TO 101
C ALLOCATE WORK AREAS FOR IN-CORE TRANSFORMATION ROUTINE TRACR:
          IF(ISR.EQ.ISS) THEN
            NBPQ=(NBP+NBP**2)/2
            NBRS=(NBR+NBR**2)/2
          ELSE
            NBPQ=NBP*NBQ
            NBRS=NBR*NBS
          ENDIF
          NAVX=NAR*NAS
          NX1MX=MAX(INTBUF,NBP*NAQ,NBQ*NAP)
          NX2MX=MAX(NBR*NBS,NAP*NAQ)
          NX3MX=MAX(NBR*NAS,NAR*NBS,NBP*NBQ)
          NVXPQ=NAVX*NBPQ
          CALL mma_allocate(X1,NX1MX,Label='X1')
          CALL mma_allocate(X2,NX2MX,Label='X2')
          CALL mma_allocate(X3,NX3MX,Label='X3')
          CALL mma_allocate(VXPQ,NVXPQ,Label='VXPQ')
          CALL TRACR( INTBUF,CMO1,CMO2,NGAM2,TUVX,
     *                X1,X2,X3,VXPQ)
          CALL mma_deallocate(X1)
          CALL mma_deallocate(X2)
          CALL mma_deallocate(X3)
          CALL mma_deallocate(VXPQ)
101      CONTINUE
102     CONTINUE
103    CONTINUE
104   CONTINUE
      Call GADSum(TUVX,NGAM2)
      END SUBROUTINE TRAINT
