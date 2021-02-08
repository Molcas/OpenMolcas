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
* Copyright (C) 1989, Per Ake Malmqvist                                *
************************************************************************
*****************************************************************
*  PROGRAM RASSI        PER-AAKE MALMQVIST
*  SUBROUTINE TRACR     IBM-3090 RELEASE 89 01 30
*  TRANSFORM ONE SPECIFIC SYMMETRY BLOCK OF TWO-ELECTRON
*  INTEGRALS. THIS ROUTINE IS CALLED FROM TRAINT.
*****************************************************************
      SUBROUTINE TRACR(LBUF,CMO1,CMO2,NGAM2,TUVX,X1,X2,X3,VXPQ)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION X1(NX1MX),X2(NX2MX),X3(NX3MX),VXPQ(NVXPQ)
      DIMENSION CMO1(NCMO),CMO2(NCMO),TUVX(NGAM2)
#include "rassi.fh"
      COMMON/TRNSFRM/ NPQ,NBPQ,NBRS,LADX,NAVX,NAP,NAQ,NAR,NAS,
     *               NBP,NBQ,NBR,NBS,ISP,ISQ,ISR,ISS,IAPR(8),
     *               LMOP1,LMOQ1,LMOR1,LMOS1,NX1MX,NX2MX,
     *               NX3MX,NVXPQ
C START LOOP OVER ORDERED AO-INTEGRALS: NPQ PQ-PAIRS IN EACH BUFFER.
C FOR EACH PQ PAIR, THERE IS A MATRIX CONTAINING THE (PQ,RS)
C INTEGRALS. INDEX S RUNS FASTEST, SO BY FORTRAN RULES, IF REGARDED
C AS A MATRIX JPQ(R,S), THEN THE MATRIX IS TRANSPOSED. IF ISR=ISS,
C THE MATRIX IS STORED IN ROW-MAJOR UNDERTRIANGULAR FORMAT.
C IN THE COMMENTS, NOTE THAT T,U, ETC DENOTES ACTUAL ORBITAL
C INDICES, WHILE IT,IU ETC ARE COUNTERS WITHIN THE SYMMETRY BLOCKS.
      IRC=0
      IOPT=1
      IPQ=0
      LPQ=0
      NPQ=0
      IRSST=1-NBRS
      DO 10 NP=1,NBP
       NQM=NBQ
       IF(ISP.EQ.ISQ) NQM=NP
       DO 9 NQ=1,NQM
        IPQ=IPQ+1
C IF NECESSARY, READ IN A FRESH INTEGRAL BUFFER OF NPQ MATRICES:
        IF(LPQ.EQ.NPQ) THEN
          CALL RDORD(IRC,IOPT,ISP,ISQ,ISR,ISS,X1,LBUF,NPQ)
          IOPT=2
          LPQ=0
          IRSST=1-NBRS
        ENDIF
        LPQ=LPQ+1
        IRSST=IRSST+NBRS
C START TRANSFORMATION OF THIS PQ=PAIR
        IF(ISR.EQ.ISS) THEN
          CALL SQUARE(X1(IRSST),X2,1,NBR,NBR)
        ELSE
          CALL DCOPY_(NBR*NBS,X1(IRSST),1,X2,1)
        ENDIF
C X2 CONTAINS THE MATRIX JPQ TRANSPOSED, X2(IS,IR)=JPQ(R,S)=(PQ/RS).
*         CALL MXMA(X2,        NBS,1,
*     *             CMO2(LMOS1),1,NBS,
*     *             X3,         1,NBR,
*     *             NBR,NBS,NAS)
*         CALL MXMA(X3,         NBR,1,
*     *             CMO1(LMOR1),1,NBR,
*     *             X2,         1,NAS,
*     *             NAS,NBR,NAR)
         CALL DGEMM_('T','N',NBR,NAS,NBS,1.0D0,
     &              X2,NBS,CMO2(LMOS1),NBS,
     &        0.0D0,X3,NBR)
         CALL DGEMM_('T','N',NAS,NAR,NBR,1.0D0,
     &              X3,NBR,CMO1(LMOR1),NBR,
     &        0.0D0,X2,NAS)
C X2 IS TRANSFORMED JPQ MATRIX, TRANSPOSED: X2(IX,IV)=(PQ/VX).
C SORT THE MATRIX X2 INTO VXPQ (SORT AFTER PQ INSTEAD OF VX).
        CALL DCOPY_(NAVX,X2,1,VXPQ(IPQ),NBPQ)
9      CONTINUE
10    CONTINUE
C FIRST HALF TRANSFORMATION IS NOW DONE.VXPQ CONTAINS HALF TRANS-
C FORMED INTEGRALS FOR THIS SYMMETRY BLOCK: VXPQ(IPQ,IX,IV)=(PQ/VX).
C NOTE: V IS FIRST ORBITAL SET, X SECOND
C NOW TRANSFORM INDICES PQ TO TU FOR ALL (TU).GE.(VX)
C FORTRAN ORDER IMPLIES THEN THAT U>=X AND, IF U=X, THEN T>=V.
      IPQST=1-NBPQ
      DO 20 IV=1,NAR
        IVF=IV+IAPR(ISR)
        DO 21 IX=1,NAS
          IXF=IX+IAPR(ISS)
          IVX=IVF+NASHT*(IXF-1)
          IPQST=IPQST+NBPQ
          IF(ISP.EQ.ISQ) THEN
            CALL SQUARE(VXPQ(IPQST),X3,1,NBP,NBP)
          ELSE
            CALL DCOPY_(NBPQ,VXPQ(IPQST),1,X3,1)
          ENDIF
C X3 IS HALF-TRANSFORMED JVX MATRIX, TRANSPOSED. X3(IQ,IP)=(PQ/VX).
C WHEN TRANSFORMING, SKIP INDICES U=1..IUM.
          IUM=0
          IF(ISQ.EQ.ISS) IUM=IX-1
*          CALL MXMA(X3,                 NBQ,1,
*     *              CMO2(LMOQ1+IUM*NBQ),1,NBQ,
*     *              X1,                 1,NBP,
*     *              NBP,NBQ,NAQ-IUM)
         CALL DGEMM_('T','N',NBP,NAQ-IUM,NBQ,1.0D0,
     &              X3,NBQ,CMO2(LMOQ1+IUM*NBQ),NBQ,
     &        0.0D0,X1,NBP)
C X1(IU,IP) = (PU/VX) FOR THE GIVEN VX, ALL P, AND U>=X, WHERE P
C IS BASIS FUNCTIONS IN SYMMETRY ISP, U AND X ARE ACTIVE ORBITALS
C OF STATE 2 IN SYMMETRIES ISQ AND ISS, RESP., AND V IS ACTIVE ORB
C OF STATE 1 AND HAS SYMMETRY ISR.
*          CALL MXMA(X1,         NBP,1,
*     *              CMO1(LMOP1),1,NBP,
*     *              X2,         1,NAQ-IUM,
*     *              NAQ-IUM,NBP,NAP)
         CALL DGEMM_('T','N',NAQ-IUM,NAP,NBP,1.0D0,
     &              X1,NBP,CMO1(LMOP1),NBP,
     &        0.0D0,X2,NAQ-IUM)
C X2 NOW CONTAINS THE TRANSFORMED JVX MATRIX TRANSPOSED, I.E.,
C X2(IU-IUM,IT)=JVX(T,U)=(TU/VX), IU=IUM+1,..,NAQ,  IT=1,..,NAP.
          II=0
          DO 14 IT=1,NAP
            ITF=IT+IAPR(ISP)
            DO 15 IU=IUM+1,NAQ
              IUF=IU+IAPR(ISQ)
              II=II+1
              ITU=ITF+NASHT*(IUF-1)
              IF(ITU.LT.IVX) THEN
                ITUVX=(IVX*(IVX-1))/2+ITU
              ELSE
                ITUVX=(ITU*(ITU-1))/2+IVX
              END IF
              TUVX(ITUVX)=X2(II)
15          CONTINUE
14        CONTINUE
      IF(ISP.EQ.ISQ) GOTO 21
C X3 STILL CONTAINS THE JVX MATRIX, TRANSPOSED. X3(IQ,IP)=(PQ/VX).
C PERM SYMMETRY OF P AND Q MEANS WE CAN LET THEM SWITCH IDENTITY.
C FOR REMAINDER OF LOOP 21, IF ISP.NE.ISQ, THEN REGARD X3 AS
C CONTAINING X3(IQ,IP)=(QP/VX) AND TRANSFORM AS BEFORE.
C WHEN TRANSFORMING, SKIP INDICES U=1..IUM.
          IUM=0
          IF(ISP.EQ.ISS) IUM=IX-1
*          CALL MXMA(X3,                 1,NBQ,
*     *              CMO2(LMOP1+IUM*NBP),1,NBP,
*     *              X1,                 1,NBQ,
*     *              NBQ,NBP,NAP-IUM)
         CALL DGEMM_('N','N',NBQ,NAP-IUM,NBP,1.0D0,
     &              X3,NBQ,CMO2(LMOP1+IUM*NBP),NBP,
     &        0.0D0,X1,NBQ)
C X1(IQ,IU) = (QU/VX) FOR THE GIVEN VX, ALL Q, AND U>=X, WHERE Q
C IS BASIS FUNCTIONS IN SYMMETRY ISQ, U AND X ARE ACTIVE ORBITALS
C OF STATE 2 IN SYMMETRIES ISP AND ISS, RESP., AND V IS ACTIVE ORB
C OF STATE 1 AND HAS SYMMETRY ISR.
*          CALL MXMA(X1,         NBQ,1,
*     *              CMO1(LMOQ1),1,NBQ,
*     *              X2,         1,NAP-IUM,
*     *              NAP-IUM,NBQ,NAQ)
         CALL DGEMM_('T','N',NAP-IUM,NAQ,NBQ,1.0D0,
     &              X1,NBQ,CMO1(LMOQ1),NBQ,
     &        0.0D0,X2,NAP-IUM)
C X2 NOW CONTAINS THE TRANSFORMED JVX MATRIX TRANSPOSED, I.E.,
C X2(IU-IUM,IT)=JVX(T,U)=(TU/VX), IU=IUM+1,..,NAP,  IT=1,..,NAQ.
          II=0
          DO 16 IT=1,NAQ
            ITF=IT+IAPR(ISQ)
            DO 17 IU=IUM+1,NAP
              IUF=IU+IAPR(ISP)
              II=II+1
              ITU=ITF+NASHT*(IUF-1)
              IF(ITU.LT.IVX) THEN
                ITUVX=(IVX*(IVX-1))/2+ITU
              ELSE
                ITUVX=(ITU*(ITU-1))/2+IVX
              END IF
              TUVX(ITUVX)=X2(II)
17          CONTINUE
16        CONTINUE
21      CONTINUE
20    CONTINUE
      IF(ISR.EQ.ISS) GOTO 100
C NOW REPEAT IT ALL OVER AGAIN. THIS TIME, USE PERMUTATION
C SYMMETRY TO SWITCH IDENTITY OF BASIS INDICES R AND S.
      IRC=0
      IOPT=1
      IPQ=0
      LPQ=0
      NPQ=0
      IRSST=1-NBRS
      DO 30 NP=1,NBP
       NQM=NBQ
       IF(ISP.EQ.ISQ) NQM=NP
       DO 29 NQ=1,NQM
        IPQ=IPQ+1
C IF NECESSARY, READ IN A FRESH INTEGRAL BUFFER OF NPQ MATRICES:
        IF(LPQ.EQ.NPQ) THEN
          CALL RDORD(IRC,IOPT,ISP,ISQ,ISR,ISS,X1,LBUF,NPQ)
          IOPT=2
          LPQ=0
          IRSST=1-NBRS
        ENDIF
        LPQ=LPQ+1
        IRSST=IRSST+NBRS
        CALL DCOPY_(NBR*NBS,X1(IRSST),1,X2,1)
C X2 CONTAINS THE MATRIX JPQ , X2(IS,IR)=JPQ(S,R)=(PQ/SR).
C NOTE THAT SYMMETRY OF R IS ISS AND SYMMETRY OF S IS ISR NOW.
*         CALL MXMA(X2,         1,NBS,
*     *             CMO2(LMOR1),1,NBR,
*     *             X3,         1,NBS,
*     *             NBS,NBR,NAR)
*         CALL MXMA(X3,         NBS,1,
*     *             CMO1(LMOS1),1,NBS,
*     *             X2,         1,NAR,
*     *             NAR,NBS,NAS)
         CALL DGEMM_('N','N',NBS,NAR,NBR,1.0D0,
     &              X2,NBS,CMO2(LMOR1),NBR,
     &        0.0D0,X3,NBS)
         CALL DGEMM_('T','N',NAR,NAS,NBS,1.0D0,
     &              X3,NBS,CMO1(LMOS1),NBS,
     &        0.0D0,X2,NAR)
C X2 IS TRANSFORMED JPQ MATRIX, TRANSPOSED: X2(IX,IV)=(PQ/VX).
C SORT THE MATRIX X2 INTO VXPQ (SORT AFTER PQ INSTEAD OF VX).
        CALL DCOPY_(NAVX,X2,1,VXPQ(IPQ),NBPQ)
29     CONTINUE
30    CONTINUE
C AS BEFORE, EXCEPT THAT SYMMETRY OF V IS ISS, SYMMETRY OF X IS ISR.
C NOW TRANSFORM INDICES PQ TO TU FOR ALL (TU).GE.(VX)
C FORTRAN ORDER IMPLIES THEN THAT U>=X AND, IF U=X, THEN T>=V.
      IPQST=1-NBPQ
      DO 40 IV=1,NAS
        IVF=IV+IAPR(ISS)
        DO 41 IX=1,NAR
          IXF=IX+IAPR(ISR)
          IVX=IVF+NASHT*(IXF-1)
          IPQST=IPQST+NBPQ
          IF(ISP.EQ.ISQ) THEN
            CALL SQUARE(VXPQ(IPQST),X3,1,NBP,NBP)
          ELSE
            CALL DCOPY_(NBPQ,VXPQ(IPQST),1,X3,1)
          ENDIF
C X3 IS HALF-TRANSFORMED JVX MATRIX, TRANSPOSED. X3(IQ,IP)=(PQ/VX).
C WHEN TRANSFORMING, SKIP INDICES U=1..IUM.
          IUM=0
          IF(ISQ.EQ.ISR) IUM=IX-1
*          CALL MXMA(X3,                 NBQ,1,
*     *              CMO2(LMOQ1+IUM*NBQ),1,NBQ,
*     *              X1,                 1,NBP,
*     *              NBP,NBQ,NAQ-IUM)
          CALL DGEMM_('T','N',NBP,NAQ-IUM,NBQ,1.0D0,
     &               X3,NBQ,CMO2(LMOQ1+IUM*NBQ),NBQ,
     &        0.0D0, X1,NBP)
C X1(IU,IP) = (PU/VX) FOR THE GIVEN VX, ALL P, AND U>=X, WHERE P
C IS BASIS FUNCTIONS IN SYMMETRY ISP, U AND X ARE ACTIVE ORBITALS
C OF STATE 2 IN SYMMETRIES ISQ AND ISR, RESP., AND V IS ACTIVE ORB
C OF STATE 1 AND HAS SYMMETRY ISS.
*          CALL MXMA(X1,         NBP,1,
*     *              CMO1(LMOP1),1,NBP,
*     *              X2,         1,NAQ-IUM,
*     *              NAQ-IUM,NBP,NAP)
          CALL DGEMM_('T','N',NAQ-IUM,NAP,NBP,1.0D0,
     &               X1,NBP,CMO1(LMOP1),NBP,
     &        0.0D0, X2,NAQ-IUM)
C X2 NOW CONTAINS THE TRANSFORMED JVX MATRIX TRANSPOSED, I.E.,
C X2(IU,IT)=JVX(T,U)=(TU/VX), IU=IUM+1,..,NAQ,  IT=1,..,NAP.
          II=0
          DO 34 IT=1,NAP
            ITF=IT+IAPR(ISP)
            DO 35 IU=IUM+1,NAQ
              IUF=IU+IAPR(ISQ)
              II=II+1
              ITU=ITF+NASHT*(IUF-1)
              IF(ITU.LT.IVX) THEN
                ITUVX=(IVX*(IVX-1))/2+ITU
              ELSE
                ITUVX=(ITU*(ITU-1))/2+IVX
              END IF
              TUVX(ITUVX)=X2(II)
35          CONTINUE
34        CONTINUE
      IF(ISP.EQ.ISQ) GOTO 41
C X3 STILL CONTAINS THE JVX MATRIX, TRANSPOSED. X3(IQ,IP)=(PQ/VX).
C PERM SYMMETRY OF P AND Q MEANS WE CAN LET THEM SWITCH IDENTITY.
C FOR REMAINDER OF LOOP 41, IF ISP.NE.ISQ, THEN REGARD X3 AS
C CONTAINING X3(IQ,IP)=(QP/VX) AND TRANSFORM AS BEFORE.
C WHEN TRANSFORMING, SKIP INDICES U=1..IUM.
          IUM=0
          IF(ISP.EQ.ISR) IUM=IX-1
*          CALL MXMA(X3,                 1,NBQ,
*     *              CMO2(LMOP1+IUM*NBP),1,NBP,
*     *              X1,                 1,NBQ,
*     *              NBQ,NBP,NAP-IUM)
          CALL DGEMM_('N','N',NBQ,NAP-IUM,NBP,1.0D0,
     &               X3,NBQ,CMO2(LMOP1+IUM*NBP),NBP,
     &        0.0D0, X1,NBQ)
C X1(IQ,IU) = (QU/VX) FOR THE GIVEN VX, ALL Q, AND U>=X, WHERE Q
C IS BASIS FUNCTIONS IN SYMMETRY ISQ, U AND X ARE ACTIVE ORBITALS
C OF STATE 2 IN SYMMETRIES ISP AND ISR, RESP., AND V IS ACTIVE ORB
C OF STATE 1 AND HAS SYMMETRY ISS.
*          CALL MXMA(X1,         NBQ,1,
*     *              CMO1(LMOQ1),1,NBQ,
*     *              X2,         1,NAP-IUM,
*     *              NAP-IUM,NBQ,NAQ)
          CALL DGEMM_('T','N',NAP-IUM,NAQ,NBQ,1.0D0,
     &               X1,NBQ,CMO1(LMOQ1),NBQ,
     &        0.0D0, X2,NAP-IUM)
C X2 NOW CONTAINS THE TRANSFORMED JVX MATRIX TRANSPOSED, I.E.,
C X2(IU,IT)=JVX(T,U)=(TU/VX), IU=IUM+1,..,NAP,  IT=1,..,NAQ.
          II=0
          DO 36 IT=1,NAQ
            ITF=IT+IAPR(ISQ)
            DO 37 IU=IUM+1,NAP
              IUF=IU+IAPR(ISP)
              II=II+1
              ITU=ITF+NASHT*(IUF-1)
              IF(ITU.LT.IVX) THEN
                ITUVX=(IVX*(IVX-1))/2+ITU
              ELSE
                ITUVX=(ITU*(ITU-1))/2+IVX
              END IF
              TUVX(ITUVX)=X2(II)
37          CONTINUE
36        CONTINUE
41      CONTINUE
40    CONTINUE
100   CONTINUE
*
      RETURN
      END
