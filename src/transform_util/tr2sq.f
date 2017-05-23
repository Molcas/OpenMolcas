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
* Copyright (C) 1987, Bjorn O. Roos                                    *
*               1992, Per Ake Malmqvist                                *
************************************************************************
*--------------------------------------------*
* 1987  B. O. ROOS                           *
* DEPARTMENT OF THEORETICAL CHEMISTRY        *
* UNIVERSITY OF LUND, SWEDEN                 *
*--------------------------------------------*
      SUBROUTINE TR2Sq(CMO,X1,X2,X3,URPQ,RUPQ,TUPQ,lBuf)
C
C SECOND ORDER TWO-ELECTRON TRANSFORMATION ROUTINE
C
C THIS ROUTINE IS CALLED FOR EACH SYMMETRY BLOCK OF INTEGRALS
C (ISP,ISQ,ISR,ISS) WITH ISP.GE.ISQ AND ISR.GE.ISS.
C P,Q,R,S are SO indices.
C A,B are MO indices, counting only non-frozen and non-deleted.
C T,U are occupied MO indices, only non-frozen and non-deleted.
C INTEGRALS (AB/TU) ARE ALWAYS GENERATED
C EXCHANGE INTEGRALS (AT/BU) ARE GENERATED AS FOLLOWS:
C (AT/BU) IF ISP.GE.ISR
C (AT/UB) IF ISP.GT.ISS AND ISP.NE.ISQ
C (TA/BU) IF ISQ.GT.ISR AND ISP.NE.ISQ
C (TA/UB) IF ISQ.GE.ISS AND ISP.NE.ISQ
C
C     ********** IBM-3090 RELEASE 87 09 14 **********
C     Replace MXMA with DGEMM. P-AA Malmqvist 1992-05-06.
C
      IMPLICIT REAL*8 (A-H,O-Z)
CPAM98      COMMON/INTTRA/ISP,ISQ,ISR,ISS,NBP,NBQ,NBR,NBS,NBPQ,NBRS,IRRST,
CPAM98     &              NOCP,NOCQ,NOCR,NOCS,NPQ,LADX,LRUPQ,LURPQ,LTUPQ,
CPAM98     &              NOP,NOQ,NOR,NOS,LMOP,LMOQ,LMOR,LMOS,LMOP2,LMOQ2,
CPAM98     &              LMOR2,LMOS2,IAD13,ITP,ITQ,ITR,ITS

#include "rasdim.fh"
#include "caspt2.fh"
#include "intgrl.fh"

#include "SysDef.fh"
#include "trafo.fh"
      DIMENSION CMO(NCMO)
      DIMENSION X1(*),X2(*),X3(*),RUPQ(*),URPQ(*),TUPQ(*)

      CALL QENTER('Tr2Sq')

      Kount=0
      NSYMP=(NSYM**2+NSYM)/2
      NORU=NBR*NOCS
      NOUR=NBS*NOCR
      NOTU=NOCR*NOCS
      IF(ISR.EQ.ISS) NOTU=(NOCR**2+NOCR)/2
C
C     CHECK FOR IN CORE OR OUT OF CORE TRANSFORMATION
C
C     1. SORT OF PARTIALLY TRANSFORMED INTEGRALS (RU/PQ) ON UNIT LUHLF1
      IPQMX1=NBPQ
      IF(NBPQ*NORU.GT.LRUPQ) THEN
       IPQMX1=LRUPQ/NORU
c      WRITE(*,*) 'OUT OF CORE SORT FOR INTEGRALS (RU/PQ)',IPQMX1
       IAD1S=0
       CALL dDAFILE(LUHLF1,0,RUPQ,IPQMX1,IAD1S)
      ENDIF
      IAD1=0
      IOUT1=0
C     2. SORT OF PARTIALLY TRANSFORMED INTEGRALS (UR/PQ) ON UNIT LUHLF2
      IPQMX2=NBPQ
      IF(NBPQ*NOUR.GT.LURPQ) THEN
       IPQMX2=LURPQ/NOUR
c      WRITE(*,*) 'OUT OF CORE SORT FOR INTEGRALS (UR/PQ)',IPQMX2
       IAD2S=0
       CALL dDAFILE(LUHLF2,0,URPQ,IPQMX2,IAD2S)
      ENDIF
      IAD2=0
      IOUT2=0
C     3. SORT OF PARTIALLY TRANSFORMED INTEGRALS (TU/PQ) ON UNIT LUHLF3
      IPQMX3=NBPQ
      IF(NBPQ*NOTU.GT.LTUPQ) THEN
       IPQMX3=LTUPQ/NOTU
c      WRITE(*,*) 'OUT OF CORE SORT FOR INTEGRALS (TU/PQ)',IPQMX3
       IAD3S=0
       CALL dDAFILE(LUHLF3,0,TUPQ,IPQMX3,IAD3S)
      ENDIF
      IAD3=0
      IOUT3=0
C
C     START LOOP OVER SORTED AO-INTEGRALS: NPQ PQ-PAIRS IN EACH BUFFER
C
      IPQ=0
      LPQ=0
      NPQ=0
      iOpt=1
      iRc=0
      IRSST=1-NBRS
      DO 11 NP=1,NBP
       NQM=NBQ
       IF(ISP.EQ.ISQ) NQM=NP
       DO 10 NQ=1,NQM
        IPQ=IPQ+1
        IOUT1=IOUT1+1
        IOUT2=IOUT2+1
        IOUT3=IOUT3+1
C
C        READ IN A BLOCK OF INTEGRALS FOR NPQ PQ-VALUES
C
        IF(LPQ.EQ.NPQ) THEN
         Call RdOrd(iRc,iOpt,isP,isQ,isR,isS,X1,lBuf,nPQ)
         IF(IRC.GT.1) THEN
           WRITE(6,*)' ERROR RETURN CODE IRC=',IRC
           WRITE(6,*)' FROM RDORD, CALLED FROM TRA2.'
           CALL Abend
         END IF
         iOpt=2
         LPQ=0
         IRSST=1-NBRS
        ENDIF
        LPQ=LPQ+1
        IRSST=IRSST+NBRS
C
C       START TRANFORMATION OF THIS PQ PAIR
C
        IF(ISR.EQ.ISS) THEN
         CALL SQUARE(X1(IRSST),X2,1,NBS,NBS)
        ELSE
         CALL DCOPY_(NBRS,X1(IRSST),1,X2,1)
        ENDIF
C
C       INTEGRALS (PQ/UR)
C
        IF((ISP.NE.ISQ.AND.ISP.GT.ISS).AND.NOCR.NE.0) THEN
         CALL DGEMM_('N','N',
     &               NBS,NOCR,NBR,
     &               1.0d0,X2,NBS,
     &               CMO(LMOR2),NBR,
     &               0.0d0,X3,NBS)
C
C        SORT THESE INTEGRALS AS (UR/PQ)
C
         IF(IOUT2.GT.IPQMX2) THEN
          IOUT2=1
cvv          DO 2 I=1,NOUR
cvv           CALL dDAFILE(LUHLF2,1,URPQ(1+IPQMX2*(I-1)),IPQMX2,IAD2)
cvv    2     CONTINUE
         CALL dDAFILE(LUHLF2,1,URPQ,IPQMX2*NOUR,IAD2)
         ENDIF
         CALL DCOPY_(NOUR,X3,1,URPQ(IOUT2),IPQMX2)
        ENDIF
C
C       INTEGRALS (PQ/RU)
C
        IF(NOCS.NE.0) THEN
         CALL DGEMM_('T','N',
     &               NBR,NOCS,NBS,
     &               1.0d0,X2,NBS,
     &               CMO(LMOS2),NBS,
     &               0.0d0,X3,NBR)
C
C        SORT THESE INTEGRALS AS (RU/PQ)
C
         IF(ISP.GE.ISR) THEN
          IF(IOUT1.GT.IPQMX1) THEN
           IOUT1=1
cvv           DO 4 I=1,NORU
cvv            CALL dDAFILE(LUHLF1,1,RUPQ(1+IPQMX1*(I-1)),IPQMX1,IAD1)
cvv    4      CONTINUE
        CALL dDAFILE(LUHLF1,1,RUPQ,IPQMX1*NORU,IAD1)
          ENDIF
          CALL DCOPY_(NORU,X3,1,RUPQ(IOUT1),IPQMX1)
         ENDIF
        ENDIF
C
C       INTEGRALS (PQ/TU)
C
        IF(NOCR*NOCS.NE.0) THEN
         IF(ISR.EQ.ISS) THEN
          CALL MXMT(X3,        NBR,1,
     &              CMO(LMOR2),1,NBR,
     &              X2,
     &              NOCR,NBR)
         ELSE
          CALL DGEMM_('T','N',
     &                NOCS,NOCR,NBR,
     &                1.0d0,X3,NBR,
     &                CMO(LMOR2),NBR,
     &                0.0d0,X2,NOCS)
         ENDIF
C
C        SORT INTEGRALS (PQ/TU) INTO TUPQ (SORT AFTER PQ INSTEAD OF TU)
C
         IF(IOUT3.GT.IPQMX3) THEN
          IOUT3=1
cvv          DO 6 I=1,NOTU
cvv           CALL dDAFILE(LUHLF3,1,TUPQ(1+IPQMX3*(I-1)),IPQMX3,IAD3)
cvv    6     CONTINUE
        CALL dDAFILE(LUHLF3,1,TUPQ,IPQMX3*NOTU,IAD3)

         ENDIF
         CALL DCOPY_(NOTU,X2,1,TUPQ(IOUT3),IPQMX3)
        ENDIF
C
C       WE NOW HAVE THREE SETS OF PARTIALLY TRANSFORMED INTEGRALS
C       IN TUPQ: (TU/PQ)   TRIANGULAR FOR ISR.EQ.ISS
C       IN RUPQ: (RU/PQ)   IF ISP.GE.ISR
C       IN URPQ: (UR/PQ)   IF ISP.GT.ISQ.AND.ISP.GT.ISS
C
   10  CONTINUE
   11 CONTINUE
C
C     EMPTY LAST BUFFERS
C
      IF(IPQMX1.LT.NBPQ) THEN
cvv       DO 12 I=1,NORU
cvv        CALL dDAFILE(LUHLF1,1,RUPQ(1+IPQMX1*(I-1)),IPQMX1,IAD1)
cvv   12  CONTINUE
      CALL dDAFILE(LUHLF1,1,RUPQ,IPQMX1*NORU,IAD1)
      ENDIF
      IF(IPQMX2.LT.NBPQ) THEN
cvv       DO 13 I=1,NOUR
cvv        CALL dDAFILE(LUHLF2,1,URPQ(1+IPQMX2*(I-1)),IPQMX2,IAD2)
cvv   13  CONTINUE
      CALL dDAFILE(LUHLF2,1,URPQ,IPQMX2*NOUR,IAD2)
      ENDIF
      IF(IPQMX3.LT.NBPQ) THEN
cvv       DO 14 I=1,NOTU
cvv        CALL dDAFILE(LUHLF3,1,TUPQ(1+IPQMX3*(I-1)),IPQMX3,IAD3)
cvv   14  CONTINUE
      CALL dDAFILE(LUHLF3,1,TUPQ,IPQMX3*NOTU,IAD3)

      ENDIF
C
C     FIRST PARTIAL TRANSFORMATION FINISHED
C     SORTED INTEGRALS ARE ON UNITS LUHLF1 (RUPQ), LUHLF2 (URPQ),
C     AND LUHLF3 (TUPQ), CONTROLLED BY THE ADRESSES IAD1,IAD2, AND IAD3,
C     OR IN CORE (RUPQ, URPQ, AND TUPQ)
C
C     SECOND HALF TRANSFORMATION FOR INTEGRALS (PQ/TU)
C     FIRST SAVE THE START ADDRESS ON LUINTM FOR THIS BLOCK OF INTEGRALS
C     NOTE THAT THE SYMMETRY LABEL ISPQRS ASSUMES THAT SYMMETRY LOOPS
C     IN THE ORDER T,U,A,B FOR ALL INTEGRAL TYPES.
C
      IF(NOCR*NOCS.EQ.0) GO TO 21
      ISPQRS=((ISR**2-ISR)/2+ISS-1)*NSYMP+(ISP**2-ISP)/2+ISQ
      IAD2M(1,ISPQRS)=IAD13
      ITU=0
      DO 20 NT=1,NOCR
       NUM=NT
       IF(ISS.NE.ISR) NUM=NOCS
      DO 20 NU=1,NUM
       IPQST=1+NBPQ*ITU
       ITU=ITU+1
C
C      READ ONE BUFFER OF INTEGRALS BACK INTO CORE
C
       IF(IPQMX3.LT.NBPQ) THEN
**sta0830
          Call RBuf_tra2(LUHLF3,TUPQ,NBPQ,IPQMX3,NOTU,ITU,IPQST,IAD3S)
**end0830
       ENDIF
C
       IF(ISP.EQ.ISQ) THEN
        CALL SQUARE(TUPQ(IPQST),X2,1,NBQ,NBQ)
        CALL DGEMM_('N','N',
     &              NBQ,NOP,NBP,
     &              1.0d0,X2,NBQ,
     &              CMO(LMOP),NBP,
     &              0.0d0,X1,NBP)
        CALL MXMT(X1,       NBQ,1,
     &            CMO(LMOQ),1,NBQ,
     &            X2,
     &            NOP,NBQ)
        IX2=(NOP+NOP**2)/2
       ELSE
        CALL DGEMM_('N','N',
     &              NBQ,NOP,NBP,
     &              1.0d0,TUPQ(IPQST),NBQ,
     &              CMO(LMOP),NBP,
     &              0.0d0,X1,NBQ)
        CALL DGEMM_('T','N',
     &              NOQ,NOP,NBQ,
     &              1.0d0,CMO(LMOQ),NBQ,
     &              X1,NBQ,
     &              0.0d0,X2,NOQ)
        IX2=NOP*NOQ
       ENDIF
C
C      WRITE INTEGRALS (AB/TU) ON OUTPUT UNIT LUINTM
C      INTEGRALS FOR SYMMETRY BLOCK (ISP,ISQ,ISR,ISS) ARE STORED
C      ONE BLOCK FOR EACH TU STARTING AT ADDRESS IAD2M(1,ISPQRS).
C      TRIANGULAR IN AB AND TU IF ISP.EQ.ISQ ( AND ISR.EQ.ISS)
C
       Call GADSum(X2,IX2)
       CALL dDAFILE(LUINTM,1,X2,IX2,IAD13)
C
C      EXTRACT INTEGRALS WITH ALL INDICES ACTIVE INTO TUVX
C      Not used with caspt2. If wanted, TUVX must be declared.
C
C      IF(ISP.GE.ISR.AND.NASH(ISP)*NASH(ISQ).NE.0)
C    & CALL GTUVX(X2,TUVX,NT,NU,ITP,ITQ,ITR,ITS,ISP,ISQ)

   20 CONTINUE
   21 CONTINUE
C
C     SECOND PARTIAL TRANSFORMATION FOR INTEGRALS (PQ/RU)-> (AT/BU)
C     IF ISP.EQ.ISR THEN T.GE.U BUT ALWAYS ALL A AND B
C
      NOTU=NOCQ*NOCS
      IF(ISQ.EQ.ISS) NOTU=(NOCQ**2+NOCQ)/2
      IF(ISP.GE.ISR.AND.NOTU.NE.0) THEN
       LAR=LTUPQ/NOTU
       LR=LAR/NOP
       IF(LR.GT.NBR) LR=NBR
       LAR=NOP*LR
       IAD3S=0
       CALL dDAFILE(LUHLF3,0,TUPQ,LAR,IAD3S)
       IAD3=0
       IR=0
       DO 30 NR=1,NBR
        IR=IR+1
       DO 30 NU=1,NOCS
        IRU=NBR*(NU-1)+NR
        IPQST=1+NBPQ*(IRU-1)
C
C       READ ONE BUFFER OF INTEGRALS BACK INTO CORE
C
        IF(IPQMX1.LT.NBPQ) THEN
**sta0830
          Call RBuf_tra2(LUHLF1,RUPQ,NBPQ,IPQMX1,NORU,IRU,IPQST,IAD1S)
**end0830
        ENDIF
        IF(ISP.EQ.ISQ) THEN
         CALL SQUARE(RUPQ(IPQST),X2,1,NBQ,NBQ)
        ELSE
         CALL DCOPY_(NBPQ,RUPQ(IPQST),1,X2,1)
        ENDIF
        IF(ISQ.EQ.ISS) THEN
         CALL DGEMM_('T','N',
     &               NBP,NOCQ-NU+1,NBQ,
     &               1.0d0,X2,NBQ,
     &               CMO(LMOQ2+NBQ*(NU-1)),NBQ,
     &               0.0d0,X1,NBP)
         CALL DGEMM_('T','N',
     &               NOCQ-NU+1,NOP,NBP,
     &               1.0d0,X1,NBP,
     &               CMO(LMOP),NBP,
     &               0.0d0,X2,NOCQ-NU+1)
        ELSE
         CALL DGEMM_('T','N',
     &               NBP,NOCQ,NBQ,
     &               1.0d0,X2,NBQ,
     &               CMO(LMOQ2),NBQ,
     &               0.0d0,X1,NBP)
         CALL DGEMM_('T','N',
     &               NOCQ,NOP,NBP,
     &               1.0d0,X1,NBP,
     &               CMO(LMOP),NBP,
     &               0.0d0,X2,NOCQ)
         ENDIF
C
C       INTEGRALS (AT/RU) ARE NOW IN X2 FOR ONE VALUE OF R,U AND ALL
C       VALUES OF A,T(T.GE.U IF ISQ.EQ.ISS)
C
C       SORT THESE INTEGRALS AFTER PAIR INDEX TU, IN CORE IF POSSIBLE
C       OTHERWISE USE LUHLF3 FOR TEMPORARY STORAGE (SAME POSITION AS FOR
C       INTEGRALS (TU/PQ)
C
        IF(IR.GT.LR) THEN
         IR=1
cvv         DO 24 I=1,NOTU
cvv          CALL dDAFILE(LUHLF3,1,TUPQ(1+LAR*(I-1)),LAR,IAD3)
cvv   24    CONTINUE
        CALL dDAFILE(LUHLF3,1,TUPQ,LAR*NOTU,IAD3)
        ENDIF
        NAT=0
        DO 26 NA=1,NOP
         NTM=1
         IF(ISQ.EQ.ISS) NTM=NU
        DO 26 NT=NTM,NOCQ
         ITU=NOCS*(NT-1)+NU-1
         IF(ISQ.LT.ISS) ITU=NOCQ*(NU-1)+NT-1
         IF(ISQ.EQ.ISS) ITU=(NT**2-NT)/2+NU-1
         NAT=NAT+1
         TUPQ(LAR*ITU+NOP*(IR-1)+NA)=X2(NAT)
   26   CONTINUE
   30  CONTINUE
C
C      EMPTY LAST BUFFER IF LR.LT.NBR
C
       IF(LR.LT.NBR) THEN
cvv        DO 31 I=1,NOTU
cvv         CALL dDAFILE(LUHLF3,1,TUPQ(1+LAR*(I-1)),LAR,IAD3)
cvv   31   CONTINUE
       CALL dDAFILE(LUHLF3,1,TUPQ,LAR*NOTU,IAD3)
       ENDIF
C
C      NOW TRANSFORM INDEX R TO B
C
C      FIRST COMPUTE AND SAVE START ADDRESS FOR THIS SYMMETRY BLOCK
C
       IF(ISQ.GE.ISS) THEN
        ISPQRS=((ISQ**2-ISQ)/2+ISS-1)*NSYMP+(ISP**2-ISP)/2+ISR
        IAD2M(2,ISPQRS)=IAD13
        NTMAX=NOCQ
        NUMAX=NOCS
       ELSE
        ISPQRS=((ISS**2-ISS)/2+ISQ-1)*NSYMP+(ISP**2-ISP)/2+ISR
        IAD2M(3,ISPQRS)=IAD13
        NTMAX=NOCS
        NUMAX=NOCQ
       ENDIF
       KKTU=0
       IST=1-NOP*NBR
       DO 40 NT=1,NTMAX
        NUM=NUMAX
        IF(ISQ.EQ.ISS) NUM=NT
       DO 40 NU=1,NUM
        KKTU=KKTU+1
        IST=IST+NOP*NBR
        IF(LR.LT.NBR)THEN
**sta0830
          Call RBuf_tra2(LUHLF3,TUPQ,NBR*NOP,LAR,NOTU,KKTU,IST,IAD3S)
**end0830
        ENDIF
        CALL DGEMM_('T','T',
     &              NOR,NOP,NBR,
     &              1.0d0,CMO(LMOR),NBR,
     &              TUPQ(IST),NOP,
     &              0.0d0,X2,NOR)
C
C       WRITE THESE BLOCK OF INTEGRALS ON LUINTM
C
        Call GADSum(X2,NOP*NOR)
        CALL dDAFILE(LUINTM,1,X2,NOP*NOR,IAD13)
   40  CONTINUE
      ENDIF
C
C     SECOND PARTIAL TRANSFORMATION FOR INTEGRALS (PQ/RU)-> (TA/BU)
C
      NOTU=NOCP*NOCS
      IF((ISP.NE.ISQ.AND.ISQ.GT.ISR).AND.NOTU.NE.0) THEN
       LAR=LTUPQ/NOTU
       LR=LAR/NOQ
       IF(LR.GT.NBR) LR=NBR
       LAR=NOQ*LR
       IAD3S=0
       CALL dDAFILE(LUHLF3,0,TUPQ,LAR,IAD3S)
       IAD3=0
       IRU=0
       IR=0
       DO 50 NR=1,NBR
        IR=IR+1
       DO 50 NU=1,NOCS
        IRU=NBR*(NU-1)+NR
        IPQST=1+NBPQ*(IRU-1)
C
C       READ ONE BUFFER OF INTEGRALS BACK INTO CORE
C
        IF(IPQMX1.LT.NBPQ) THEN
**sta0830
          Call RBuf_tra2(LUHLF1,RUPQ,NBPQ,IPQMX1,NORU,IRU,IPQST,IAD1S)
**end0830
        ENDIF
        CALL DGEMM_('N','N',
     &              NBQ,NOCP,NBP,
     &              1.0d0,RUPQ(IPQST),NBQ,
     &              CMO(LMOP2),NBP,
     &              0.0d0,X1,NBQ)
        CALL DGEMM_('T','N',
     &              NOCP,NOQ,NBQ,
     &              1.0d0,X1,NBQ,
     &              CMO(LMOQ),NBQ,
     &              0.0d0,X2,NOCP)
C
C       INTEGRALS (TA/RU) ARE NOW IN X2 FOR ONE VALUE OF R,U AND ALL
C       VALUES OF T,A . NOTE THAT T AND U HERE ALWAYS ARE OF DIFFERENT
C       SYMMETRIES
C
C       SORT THESE INTEGRALS AFTER PAIR INDEX TU, IN CORE IF POSSIBLE
C       OTHERWISE USE LUHLF3 FOR TEMPORARY STORAGE (SAME POSITION AS FOR
C       INTEGRALS (TU/PQ)
C
        IF(IR.GT.LR) THEN
         IR=1
cvv         DO 44 I=1,NOTU
cvv          CALL dDAFILE(LUHLF3,1,TUPQ(1+LAR*(I-1)),LAR,IAD3)
cvv   44    CONTINUE
        CALL dDAFILE(LUHLF3,1,TUPQ,LAR*NOTU,IAD3)

        ENDIF
        NAT=0
        DO 46 NA=1,NOQ
        DO 46 NT=1,NOCP
         ITU=NOCS*(NT-1)+NU-1
         NAT=NAT+1
         TUPQ(LAR*ITU+NOQ*(IR-1)+NA)=X2(NAT)
   46   CONTINUE
   50  CONTINUE
C
C      EMPTY LAST BUFFER IF LR.LT.NBR
C
       IF(LR.LT.NBR) THEN
cvv        DO 51 I=1,NOTU
cvv         CALL dDAFILE(LUHLF3,1,TUPQ(1+LAR*(I-1)),LAR,IAD3)
cvv   51   CONTINUE
       CALL dDAFILE(LUHLF3,1,TUPQ,LAR*NOTU,IAD3)
       ENDIF
C
C      NOW TRANSFORM INDEX R TO B
C
C      FIRST COMPUTE AND SAVE START ADDRESS FOR THIS SYMMETRY BLOCK
C
       ISPQRS=((ISP**2-ISP)/2+ISS-1)*NSYMP+(ISQ**2-ISQ)/2+ISR
       IAD2M(2,ISPQRS)=IAD13
       KKTU=0
       IST=1-NOQ*NBR
       DO 60 NT=1,NOCP
       DO 60 NU=1,NOCS
        KKTU=KKTU+1
        IST=IST+NOQ*NBR
        IF(LR.LT.NBR)THEN
**sta0830
          Call RBuf_tra2(LUHLF3,TUPQ,NBR*NOQ,LAR,NOTU,KKTU,IST,IAD3S)
**end0830
        ENDIF
        CALL DGEMM_('T','T',
     &              NOR,NOQ,NBR,
     &              1.0d0,CMO(LMOR),NBR,
     &              TUPQ(IST),NOQ,
     &              0.0d0,X2,NOR)
C
C       WRITE THESE BLOCK OF INTEGRALS ON LUINTM
C
        Call GADSum(X2,NOR*NOQ)
        CALL dDAFILE(LUINTM,1,X2,NOR*NOQ,IAD13)
   60  CONTINUE
      ENDIF
C
C     SECOND PARTIAL TRANSFORMATION FOR INTEGRALS (PQ/UR)-> (AT/UB)
C
      NOTU=NOCQ*NOCR
      IF((ISP.NE.ISQ.AND.ISP.GT.ISS).AND.NOTU.NE.0) THEN
       LAR=(LRUPQ+LTUPQ)/NOTU
       LR=LAR/NOP
       IF(LR.GT.NBS) LR=NBS
       LAR=NOP*LR
       IAD3S=0
       CALL dDAFILE(LUHLF3,0,RUPQ,LAR,IAD3S)
       IAD3=0
       IRU=0
       IR=0
       DO 70 NR=1,NBS
        IR=IR+1
       DO 70 NU=1,NOCR
        IRU=NBS*(NU-1)+NR
        IPQST=1+NBPQ*(IRU-1)
C
C       READ ONE BUFFER OF INTEGRALS BACK INTO CORE
C
        IF(IPQMX2.LT.NBPQ) THEN
**sta0830
          Call RBuf_tra2(LUHLF2,URPQ,NBPQ,IPQMX2,NOUR,IRU,IPQST,IAD2S)
**end0830
        ENDIF
        CALL DGEMM_('T','N',
     &              NBP,NOCQ,NBQ,
     &              1.0d0,URPQ(IPQST),NBQ,
     &              CMO(LMOQ2),NBQ,
     &              0.0d0,X1,NBP)
        CALL DGEMM_('T','N',
     &              NOCQ,NOP,NBP,
     &              1.0d0,X1,NBP,
     &              CMO(LMOP),NBP,
     &              0.0d0,X2,NOCQ)
C
C       INTEGRALS (AT/UR) ARE NOW IN X2 FOR ONE VALUE OF R,U AND ALL
C       VALUES OF A,T. NOTE THAT T AND U HAVE DIFFERENT SYMMETRIES.
C
C       SORT THESE INTEGRALS AFTER PAIR INDEX TU, IN CORE IF POSSIBLE
C       OTHERWISE USE LUHLF3 FOR TEMPORARY STORAGE (SAME POSITION AS FOR
C       INTEGRALS (TU/PQ)
C
        IF(IR.GT.LR) THEN
         IR=1
cvv         DO 64 I=1,NOTU
cvv          CALL dDAFILE(LUHLF3,1,RUPQ(1+LAR*(I-1)),LAR,IAD3)
cvv   64    CONTINUE
        CALL dDAFILE(LUHLF3,1,RUPQ,LAR*NOTU,IAD3)
        ENDIF
        NAT=0
        DO 66 NA=1,NOP
        DO 66 NT=1,NOCQ
         ITU=NOCR*(NT-1)+NU-1
         IF(ISQ.LT.ISR) ITU=NOCQ*(NU-1)+NT-1
         NAT=NAT+1
         RUPQ(LAR*ITU+NOP*(IR-1)+NA)=X2(NAT)
   66   CONTINUE
   70  CONTINUE
C
C      EMPTY LAST BUFFER IF LR.LT.NBS
C
       IF(LR.LT.NBS) THEN
cvv        DO 71 I=1,NOTU
cvv         CALL dDAFILE(LUHLF3,1,RUPQ(1+LAR*(I-1)),LAR,IAD3)
cvv   71   CONTINUE
       CALL dDAFILE(LUHLF3,1,RUPQ,LAR*NOTU,IAD3)
       ENDIF
C
C      NOW TRANSFORM INDEX R TO B
C
C      FIRST COMPUTE AND SAVE START ADDRESS FOR THIS SYMMETRY BLOCK
C
       IF(ISQ.GE.ISR) THEN
        ISPQRS=((ISQ**2-ISQ)/2+ISR-1)*NSYMP+(ISP**2-ISP)/2+ISS
        IAD2M(2,ISPQRS)=IAD13
        NTMAX=NOCQ
        NUMAX=NOCR
       ELSE
        ISPQRS=((ISR**2-ISR)/2+ISQ-1)*NSYMP+(ISP**2-ISP)/2+ISS
        IAD2M(3,ISPQRS)=IAD13
        NTMAX=NOCR
        NUMAX=NOCQ
       ENDIF
       KKTU=0
       IST=1-NBS*NOP
       DO 80 NT=1,NTMAX
       DO 80 NU=1,NUMAX
        KKTU=KKTU+1
        IST=IST+NBS*NOP
        IF(LR.LT.NBS)THEN
**sta0830
          Call RBuf_tra2(LUHLF3,RUPQ,NBS*NOP,LAR,NOTU,KKTU,IST,IAD3S)
**end0830
        ENDIF
        CALL DGEMM_('T','T',
     &              NOS,NOP,NBS,
     &              1.0d0,CMO(LMOS),NBS,
     &              RUPQ(IST),NOP,
     &              0.0d0,X2,NOS)
C
C       WRITE THESE BLOCK OF INTEGRALS ON LUINTM
C
        CAll GADSum(X2,NOS*NOP)
        CALL dDAFILE(LUINTM,1,X2,NOS*NOP,IAD13)
   80  CONTINUE
      ENDIF
C
C     SECOND PARTIAL TRANSFORMATION FOR INTEGRALS (PQ/UR)-> (TA/UB)
C
      NOTU=NOCP*NOCR
      IF((ISP.NE.ISQ.AND.ISQ.GE.ISS).AND.NOTU.NE.0) THEN
       IF(ISP.EQ.ISR) NOTU=(NOCP**2+NOCP)/2
       LAR=(LRUPQ+LTUPQ)/NOTU
       LR=LAR/NOQ
       IF(LR.GT.NBS) LR=NBS
       LAR=NOQ*LR
       IAD3S=0
       CALL dDAFILE(LUHLF3,0,RUPQ,LAR,IAD3S)
       IAD3=0
       IRU=0
       IR=0
       DO 90 NR=1,NBS
        IR=IR+1
       DO 90 NU=1,NOCR
        IRU=NBS*(NU-1)+NR
        IPQST=1+NBPQ*(IRU-1)
C
C       READ ONE BUFFER OF INTEGRALS BACK INTO CORE
C
        IF(IPQMX2.LT.NBPQ) THEN
**sta0830
          Call RBuf_tra2(LUHLF2,URPQ,NBPQ,IPQMX2,NOUR,IRU,IPQST,IAD2S)
**end0830
        ENDIF
        IF(ISP.EQ.ISR) THEN
         CALL DGEMM_('N','N',
     &               NBQ,NOCP-NU+1,NBP,
     &               1.0d0,URPQ(IPQST),NBQ,
     &               CMO(LMOP2+NBP*(NU-1)),NBP,
     &               0.0d0,X1,NBQ)
         CALL DGEMM_('T','N',
     &               NOCP-NU+1,NOQ,NBQ,
     &               1.0d0,X1,NBQ,
     &               CMO(LMOQ),NBQ,
     &               0.0d0,X2,NOCP-NU+1)
        ELSE
         CALL DGEMM_('N','N',
     &               NBQ,NOCP,NBP,
     &               1.0d0,URPQ(IPQST),NBQ,
     &               CMO(LMOP2),NBP,
     &               0.0d0,X1,NBQ)
         CALL DGEMM_('T','N',
     &               NOCP,NOQ,NBQ,
     &               1.0d0,X1,NBQ,
     &               CMO(LMOQ),NBQ,
     &               0.0d0,X2,NOCP)
        ENDIF
C
C       INTEGRALS (TA/UR) ARE NOW IN X2 FOR ONE VALUE OF R,U AND ALL
C       VALUES OF A,T
C
C       SORT THESE INTEGRALS AFTER PAIR INDEX TU, IN CORE IF POSSIBLE
C       OTHERWISE USE LUHLF3 FOR TEMPORARY STORAGE (SAME POSITION AS FOR
C       INTEGRALS (TU/PQ)
C
        IF(IR.GT.LR) THEN
         IR=1
cvv         DO 84 I=1,NOTU
cvv          CALL dDAFILE(LUHLF3,1,RUPQ(1+LAR*(I-1)),LAR,IAD3)
cvv   84    CONTINUE
        CALL dDAFILE(LUHLF3,1,RUPQ,LAR*NOTU,IAD3)
        ENDIF
        NAT=0
        DO 86 NA=1,NOQ
         NTM=1
         IF(ISP.EQ.ISR) NTM=NU
        DO 86 NT=NTM,NOCP
         ITU=NOCR*(NT-1)+NU-1
         IF(ISP.LT.ISR) ITU=NOCP*(NU-1)+NT-1
         IF(ISP.EQ.ISR) ITU=(NT**2-NT)/2+NU-1
         NAT=NAT+1
         RUPQ(LAR*ITU+NOQ*(IR-1)+NA)=X2(NAT)
   86   CONTINUE
   90  CONTINUE
C
C      EMPTY LAST BUFFER IF LR.LT.NBS
C
       IF(LR.LT.NBS) THEN
cvv        DO 91 I=1,NOTU
cvv         CALL dDAFILE(LUHLF3,1,RUPQ(1+LAR*(I-1)),LAR,IAD3)
cvv   91   CONTINUE
       CALL dDAFILE(LUHLF3,1,RUPQ,LAR*NOTU,IAD3)
       ENDIF
C
C      NOW TRANSFORM INDEX R TO B
C
C       FIRST COMPUTE AND SAVE START ADDRESS FOR THIS SYMMETRY BLOCK
C
       IF(ISP.GE.ISR) THEN
        ISPQRS=((ISP**2-ISP)/2+ISR-1)*NSYMP+(ISQ**2-ISQ)/2+ISS
        IAD2M(2,ISPQRS)=IAD13
        NTMAX=NOCP
        NUMAX=NOCR
       ELSE
        ISPQRS=((ISR**2-ISR)/2+ISP-1)*NSYMP+(ISQ**2-ISQ)/2+ISS
        IAD2M(3,ISPQRS)=IAD13
        NTMAX=NOCR
        NUMAX=NOCP
       ENDIF
       KKTU=0
       IST=1-NOQ*NBS
       DO 100 NT=1,NTMAX
        NUM=NUMAX
        IF(ISP.EQ.ISR) NUM=NT
       DO 100 NU=1,NUM
        KKTU=KKTU+1
        IST=IST+NOQ*NBS
        IF(LR.LT.NBS)THEN
**sta0830
          Call RBuf_tra2(LUHLF3,RUPQ,NBS*NOQ,LAR,NOTU,KKTU,IST,IAD3S)
**end0830
        ENDIF
        CALL DGEMM_('T','T',
     &              NOS,NOQ,NBS,
     &              1.0d0,CMO(LMOS),NBS,
     &              RUPQ(IST),NOQ,
     &              0.0d0,X2,NOS)
C
C       WRITE THESE BLOCK OF INTEGRALS ON LUINTM
C
        Call GADSum(X2,NOS*NOQ)
        CALL dDAFILE(LUINTM,1,X2,NOS*NOQ,IAD13)
  100  CONTINUE
      ENDIF
C
C     END OF TRANSFORMATION FOR THIS SYMMETRY BLOCK
C
C     IAD2M CONTAINS START ADRESS FOR EACH TYPE OF INTEGRALS:
C     IAD2M(1,ISPQRS)   COULOMB INTEGRALS (AB|TU)
C     IAD2M(2,ISPQRS)   EXCHANGE INTEGRALS <AB|TU> FOR SYM T > SYM U
C     IAD2M(3,ISPQRS)     EXCHANGE INTEGRALS <AB|TU> FOR SYM T < SYM U
C     THE LAST ADRESS IS ZERO IF SYM T = SYM U
C     TO SEE HOW THE INTEGRALS ARE USED LOOK IN RDINT2
C

      CALL QEXIT('Tr2Sq')

      RETURN
      END
