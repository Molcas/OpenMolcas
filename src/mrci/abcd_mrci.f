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
      SUBROUTINE ABCD_MRCI(INTSYM,indx,ISAB,C,S,ACBDS,ACBDT,BUFIN)
      IMPLICIT REAL*8 (A-H,O-Z)
#include "SysDef.fh"
#include "mrci.fh"
      DIMENSION INTSYM(*),indx(*),ISAB(NVIRT,NVIRT),
     *          C(*),S(*),ACBDS(*),ACBDT(*),
     *          BUFIN(*)
*
cvv hand-made loop unrolling to fix a bug in GCC 3.x
      IAD16=0
      INSIN=KBUFF1
      CALL CSCALE(indx,INTSYM,C,SQ2)
      CALL CSCALE(indx,INTSYM,S,SQ2INV)
      NVT=IROW(NVIRT+1)
      NOV=(NVT-1)/IPASS+1
      IACMAX=0
c      DO 70 ISTEP=1,IPASS
       ISTEP=1
       if(IPASS.lt.1) goto 670

770   IACMIN=IACMAX+1
      IACMAX=IACMAX+NOV
      IF(IACMAX.GT.NVT)IACMAX=NVT
      IF(IACMIN.GT.IACMAX)GO TO 70
c      DO 40 ISYM=1,NSYM
        ISYM=1
        if(NSYM.lt.1) goto 640
740   IST1=IRC(3)+JJS(ISYM+9)+1
      IFIN1=IRC(3)+JJS(ISYM+10)
      INPS=IFIN1-IST1+1
      IST2=IRC(2)+JJS(ISYM)+1
      IFIN2=IRC(2)+JJS(ISYM+1)
      INPT=IFIN2-IST2+1
      ITAIL=INPS+INPT
      IF(ITAIL.EQ.0)GO TO 40
      IN1=-NVIRT
c      DO 50 NA=1,NVIRT
       NA=1
       if(NVIRT.lt.1) goto 650
750      IN1=IN1+NVIRT
c      DO 60 NC=1,NA
       NC=1
       if(NA.lt.1) goto 660

760    IAC=IROW(NA)+NC
      IF(IAC.LT.IACMIN)GO TO 60
      IF(IAC.GT.IACMAX)GO TO 60
      IF(NA.EQ.1)GO TO 60
      NSAC=MUL(NSM(LN+NA),NSM(LN+NC))
      NSACL=MUL(NSAC,LSYM)
      IF(NSACL.NE.ISYM)GO TO 60
      ISAC=ISAB(NA,NC)
      NSC=NSM(LN+NC)
      NDMAX=NVIRP(NSC)+NVIR(NSC)
      IF(NDMAX.GT.NA)NDMAX=NA
      INS=ISAB(NA,NDMAX)
C MOVE INS ITEMS FROM FILE, UNIT 16, VIA BUFFER, INTO ACBDS,
C AND THEN INTO ACBDT:
      ILOOP=0
72    INSB=INS
73    IF(INSIN.LT.KBUFF1)GO TO 75
C INSB ITEMS REMAIN TO MOVE.
C INSIN ITEMS HAVE ALREADY BEEN MOVED FROM THE BUFFER.
      CALL dDAFILE(Lu_80,2,BUFIN,KBUFF1,IAD16)
      INSIN=0
75    INB=KBUFF1-INSIN
C INB FRESH ITEMS ARE STILL REMAINING IN BUFFER.
      INUMB=MIN(INSB,INB)
C MOVE INUMB ITEMS.
      IST=INS-INSB+1
      IF(ILOOP.EQ.0)CALL DCOPY_(INUMB,BUFIN(INSIN+1),1,ACBDS(IST),1)
      IF(ILOOP.EQ.1)CALL DCOPY_(INUMB,BUFIN(INSIN+1),1,ACBDT(IST),1)
      INSIN=INSIN+INUMB
      INSB=INSB-INUMB
      IF(INSB.GT.0)GO TO 73
      ILOOP=ILOOP+1
      IF(ILOOP.EQ.1)GO TO 72
C INS ITEMS HAVE BEEN TRANSFERRED TO ACBDS AND TO ACBDT.
      IF(INPS.EQ.0)GO TO 11
      INDA=IST1
      if(IFIN1.lt.IST1) goto 610
c      DO 10 INDA=IST1,IFIN1
710   TERM=DDOT_(INS,C(indx(INDA)+1),1,ACBDS,1)
      S(indx(INDA)+ISAC)=S(indx(INDA)+ISAC)+TERM
      CALL DAXPY_(INS,C(indx(INDA)+ISAC),ACBDS,1,S(indx(INDA)+1),1)
c10    CONTINUE
       INDA=INDA+1
       if(INDA.le.IFIN1) goto 710
610    Continue
11    IF(INPT.EQ.0.OR.NA.EQ.NC)GO TO 60
       INDA=IST2
      if(IFIN2.lt.IST2) goto 630
c      DO 30 INDA=IST2,IFIN2
730   TERM=DDOT_(INS,C(indx(INDA)+1),1,ACBDT,1)
      S(indx(INDA)+ISAC)=S(indx(INDA)+ISAC)+TERM
      CALL DAXPY_(INS,C(indx(INDA)+ISAC),ACBDT,1,S(indx(INDA)+1),1)
c30    CONTINUE
       INDA=INDA+1
       if(INDA.le.IFIN2) goto 730
630   Continue
60    CONTINUE
       NC=NC+1
       if(NC.le.NA) goto 760
660    Continue
cvv end of unrolling loop
c        NC=NC+1
c        if(NC.le.NA) goto 61
c50    CONTINUE
        NA=NA+1
        if(NA.le.NVIRT) goto 750
650   Continue
40    CONTINUE
       ISYM=ISYM+1
       if(ISYM.le.NSYM) goto 740
640   Continue

70    CONTINUE
       ISTEP=ISTEP+1
       IF(ISTEP.le.IPASS) goto 770
670   Continue
      CALL CSCALE(indx,INTSYM,C,SQ2INV)
      CALL CSCALE(indx,INTSYM,S,SQ2)
      RETURN
      END
