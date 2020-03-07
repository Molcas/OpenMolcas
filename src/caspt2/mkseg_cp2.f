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
* Copyright (C) 1994, Per Ake Malmqvist                                *
************************************************************************
*--------------------------------------------*
* 1994  PER-AAKE MALMQUIST                   *
* DEPARTMENT OF THEORETICAL CHEMISTRY        *
* UNIVERSITY OF LUND                         *
* SWEDEN                                     *
*--------------------------------------------*
      SUBROUTINE MKSEG_CP2(IDRT,IDOWN,LTV,IVR,MVL,MVR,ISGMNT,VSGMNT)
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER(26) CC1,CC2,CTVPT,CBVPT,CSVC
      CHARACTER(20) FRML(7)
#include "pt2_guga.fh"
      COMMON /SEGTAB/ IC1(26),IC2(26),ITVPT(26),IBVPT(26),ISVC(26),
     &                NIVR,LIVR,NSGMNT,LSGMNT
      DIMENSION IVR(NVERT,2),ISGMNT(NVERT,26),VSGMNT(NVERT,26)
      DIMENSION IDRT(NVERT,5),IDOWN(NVERT,0:3),LTV(-1:NLEV)
      DIMENSION MVL(NMIDV,2),MVR(NMIDV,2)
C     PARAMETER (LTAB=1, NTAB=2, IATAB=3, IBTAB=4, ICTAB=5)
      PARAMETER (IATAB=3, IBTAB=4)
      PARAMETER (ZERO=0.0D0, ONE=1.0D00)
#ifdef _DEBUG_
      CHARACTER(20) TEXT
#endif

C PURPOSE: CONSTRUCT THE TABLES ISGMNT AND VSGMNT.
C ISGMNT(IVLT,ISGT) REFERS TO A SEGMENT OF THE SEGMENT TYPE
C    ISGT=1,..,26, WHOSE TOP LEFT VERTEX IS IVLT. ISGMNT GIVES
C    ZERO IF THE SEGMENT IS IMPOSSIBLE IN THE GRAPH DEFINED BY
C    THE PALDUS TABLE IDRT. ELSE IT IS THE BOTTOM LEFT VERTEX
C    NUMBER OF THE SEGMENT. THE SEGMENT VALUE IS THEN VSGMNT.

      CC1=  '01230201011230122313230123'
      CC2=  '01231323012230112302010123'
      CTVPT='00000000111112222211223333'
      CBVPT='00001122112112212233333333'
      CSVC ='11111615124721732215161111'
      FRML(1)='        1.0         '
      FRML(2)='       -1.0         '
      FRML(3)='        1/(B+1)     '
      FRML(4)='       -1/(B+1)     '
      FRML(5)='  SQRT(   B /(B+1)) '
      FRML(6)='  SQRT((B+2)/(B+1)) '
      FRML(7)='  SQRT(B(B+2))/(B+1)'
      READ(CC1,'(26I1)') IC1
      READ(CC2,'(26I1)') IC2
      READ(CTVPT,'(26I1)') ITVPT
      READ(CBVPT,'(26I1)') IBVPT
      READ(CSVC,'(26I1)') ISVC
      DO 5 IV=1,NVERT
        IVR(IV,1)=0
   5    IVR(IV,2)=0
      DO 20 LEV=1,NLEV
        IV1=LTV(LEV)
        IV2=LTV(LEV-1)-1
        DO 20 IVL=IV1,IV2
          IAL=IDRT(IVL,IATAB)
          IBL=IDRT(IVL,IBTAB)
          DO 10 IV=IVL+1,IV2
            IA=IDRT(IV,IATAB)
            IF(IA.EQ.IAL) THEN
              IB=IDRT(IV,IBTAB)
              IF(IB.EQ.(IBL-1)) IVR(IVL,1)=IV
            ELSE IF (IA.EQ.(IAL-1)) THEN
              IB=IDRT(IV,IBTAB)
              IF(IB.EQ.(IBL+1)) IVR(IVL,2)=IV
            END IF
  10      CONTINUE
  20  CONTINUE
C CONSTRUCT THE MVL AND MVR TABLES:
      DO 30 IVL=MIDV1,MIDV2
        MVLL=IVL-MIDV1+1
        MVRR=0
        IF(IVR(IVL,1).NE.0) MVRR=IVR(IVL,1)-MIDV1+1
        MVR(MVLL,1)=MVRR
        MVRR=0
        IF(IVR(IVL,2).NE.0) MVRR=IVR(IVL,2)-MIDV1+1
        MVR(MVLL,2)=MVRR
        MVL(MVLL,1)=0
        MVL(MVLL,2)=0
  30  CONTINUE
      DO 31 MV=1,NMIDV
        IF(MVR(MV,1).NE.0) MVL(MVR(MV,1),1)=MV
        IF(MVR(MV,2).NE.0) MVL(MVR(MV,2),2)=MV
  31  CONTINUE
#ifdef _DEBUG_
        WRITE(6,*)
        WRITE(6,*)' MIDVERT PAIR TABLES MVL,MVR IN MKSEG:'
        WRITE(6,*)' MVL TABLE:'
        WRITE(6,1234)(MV,MVL(MV,1),MVL(MV,2),MV=1,NMIDV)
        WRITE(6,*)' MVR TABLE:'
        WRITE(6,1234)(MV,MVR(MV,1),MVR(MV,2),MV=1,NMIDV)
 1234   FORMAT(3(3(1X,I4),4X))
        WRITE(6,*)
        WRITE(6,*)' VERTEX PAIR TABLE IVR IN MKSEG:'
        WRITE(6,1234)(IVL,IVR(IVL,1),IVR(IVL,2),IVL=1,NVERT)
#endif

C INITIALIZE SEGMENT TABLES, AND MARK VERTICES AS UNUSABLE:
      DO 40 IVLT=1,NVERT
        DO 40 ISGT=1,26
          ISGMNT(IVLT,ISGT)=0
  40      VSGMNT(IVLT,ISGT)=ZERO
      DO 100 IVLT=1,NVERT
        DO 100 ISGT=1,26
          ITT=ITVPT(ISGT)
          IVRT=IVLT
          IF((ITT.EQ.1).OR.(ITT.EQ.2)) IVRT=IVR(IVLT,ITT)
          IF(IVRT.EQ.0) GOTO 100
          IVLB=IDOWN(IVLT,IC1(ISGT))
          IF(IVLB.EQ.0) GOTO 100
          IVRB=IDOWN(IVRT,IC2(ISGT))
          IF(IVRB.EQ.0) GOTO 100
C SEGMENT IS NOW ACCEPTED AS POSSIBLE.

          ISGMNT(IVLT,ISGT)=IVLB
          IB=IDRT(IVLT,IBTAB)
          GOTO (1001,1002,1003,1004,1005,1006,1007) ISVC(ISGT)
 1001     V=ONE
          GOTO 99
 1002     V=-ONE
          GOTO 99
 1003     V=ONE/DBLE(1+IB)
          GOTO 99
 1004     V=-ONE/DBLE(1+IB)
          GOTO 99
 1005     V=SQRT(DBLE(IB)/DBLE(1+IB))
          GOTO 99
 1006     V=SQRT(DBLE(2+IB)/DBLE(1+IB))
          GOTO 99
 1007     V=SQRT(DBLE(IB*(2+IB)))/DBLE(1+IB)
  99      VSGMNT(IVLT,ISGT)=V
 100  CONTINUE

#ifdef _DEBUG_
        WRITE(6,*)' SEGMENT TABLE IN MKSEG.'
        WRITE(6,*)' VLT SGT ICL ICR VLB       SEGMENT TYPE    ',
     &            '     FORMULA'
        DO 200 IV=1,NVERT
          DO 200 ISGT=1,26
            ID=ISGMNT(IV,ISGT)
            IF(ID.EQ.0) GOTO 200
            ICL=IC1(ISGT)
            ICR=IC2(ISGT)
            IF(ISGT.LE.4) TEXT='  WALK SECTION.'
            IF((ISGT.GE.5).AND.(ISGT.LE.8)) TEXT=' TOP SEGMENT.'
            IF((ISGT.GE.9).AND.(ISGT.LE.18)) TEXT=' MID-SEGMENT.'
            IF((ISGT.GE.19).AND.(ISGT.LE.22)) TEXT=' BOTTOM SEGMENT.'
            IF(ISGT.GT.22) TEXT=' DOWN-WALK SECTION.'
            WRITE(6,2345) IV,ISGT,ICL,ICR,ID,TEXT,FRML(ISVC(ISGT))
 2345       FORMAT(1X,5I4,5X,A20,5X,A20)
  200     CONTINUE
#endif
      RETURN
      END
