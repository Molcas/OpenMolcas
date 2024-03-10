!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 1994, Per Ake Malmqvist                                *
!***********************************************************************
!--------------------------------------------*
! 1994  PER-AAKE MALMQUIST                   *
! DEPARTMENT OF THEORETICAL CHEMISTRY        *
! UNIVERSITY OF LUND                         *
! SWEDEN                                     *
!--------------------------------------------*
!#define _DEBUGPRINT_
      SUBROUTINE MKSEG(SGS,CIS,EXS)

! PURPOSE: CONSTRUCT THE TABLES ISGMNT AND VSGMNT.
! ISGMNT(IVLT,ISGT) REFERS TO A SEGMENT OF THE SEGMENT TYPE
!    ISGT=1,..,26, WHOSE TOP LEFT VERTEX IS IVLT. ISGMNT GIVES
!    ZERO IF THE SEGMENT IS IMPOSSIBLE IN THE GRAPH DEFINED BY
!    THE PALDUS TABLE IDRT. ELSE IT IS THE BOTTOM LEFT VERTEX
!    NUMBER OF THE SEGMENT. THE SEGMENT VALUE IS THEN VSGMNT.

      use Struct, only: SGStruct, CIStruct,EXStruct
      use stdalloc, only: mma_allocate
      IMPLICIT None
      Type (SGStruct) SGS
      Type (CIStruct) CIS
      Type (EXStruct) EXS

! local stuff
      CHARACTER(LEN=26) CC1,CC2,CTVPT,CBVPT,CSVC
#include "segtab.fh"
      Integer, PARAMETER :: IATAB=3, IBTAB=4
      Real*8, PARAMETER :: ZERO=0.0D0, ONE=1.0D00

#ifdef _DEBUGPRINT_
      CHARACTER(LEN=20) FRML(7)
      CHARACTER(LEN=20) TEXT
      Integer :: ICL, ICR, ID
#endif

      Integer:: IA, IAL, IB, IBL, ISGT, ITT, IV, IV1, IV2, IVL, IVLB,   &
     &          IVLT, IVRB, IVRT, LEV, MV, MVLL, MVRR
      Real*8 :: V
      Call mma_allocate(CIS%IVR,SGS%nVert,2,Label='CIS%IVR')
      Call mma_allocate(CIS%ISGM,SGS%nVert,26,Label='CIS%ISGM')
      Call mma_allocate(CIS%VSGM,SGS%nVert,26,Label='CIS%VSGM')
      Call mma_allocate(EXS%MVL,CIS%nMidV,2,Label='EXS%MVL')
      Call mma_allocate(EXS%MVR,CIS%nMidV,2,Label='EXS%MVR')

! Dereference SGS
      Associate (iDRT=>SGS%DRT, iDown=>SGS%Down, LTV=>SGS%LTV,          &
     &           MVSTA=>SGS%MVSta, MVEnd=>SGS%MVEnd,                    &
     &           nLev=>SGS%nLev, nMidV=>CIS%nMidV,                      &
     &           MVL=>EXS%MVL, MVR=>EXS%MVR, nVert=>SGS%nVert,          &
     &           IVR=>CIS%IVR, iSGMNT => CIS%ISGM,                      &
     &           VSGMNT => CIS%VSGM)

      CC1=  '01230201011230122313230123'
      CC2=  '01231323012230112302010123'
      CTVPT='00000000111112222211223333'
      CBVPT='00001122112112212233333333'
      CSVC ='11111615124721732215161111'
#ifdef _DEBUGPRINT_
      FRML(1)='        1.0         '
      FRML(2)='       -1.0         '
      FRML(3)='        1/(B+1)     '
      FRML(4)='       -1/(B+1)     '
      FRML(5)='  SQRT(   B /(B+1)) '
      FRML(6)='  SQRT((B+2)/(B+1)) '
      FRML(7)='  SQRT(B(B+2))/(B+1)'
#endif
      READ(CC1,'(26I1)') IC1
      READ(CC2,'(26I1)') IC2
      READ(CTVPT,'(26I1)') ITVPT
      READ(CBVPT,'(26I1)') IBVPT
      READ(CSVC,'(26I1)') ISVC
      DO IV=1,NVERT
        IVR(IV,1)=0
        IVR(IV,2)=0
      END DO
      DO LEV=1,NLEV
        IV1=LTV(LEV)
        IV2=LTV(LEV-1)-1
        DO IVL=IV1,IV2
          IAL=IDRT(IVL,IATAB)
          IBL=IDRT(IVL,IBTAB)
          DO IV=IVL+1,IV2
            IA=IDRT(IV,IATAB)
            IF(IA.EQ.IAL) THEN
              IB=IDRT(IV,IBTAB)
              IF(IB.EQ.(IBL-1)) IVR(IVL,1)=IV
            ELSE IF (IA.EQ.(IAL-1)) THEN
              IB=IDRT(IV,IBTAB)
              IF(IB.EQ.(IBL+1)) IVR(IVL,2)=IV
            END IF
          END DO
        END DO
      END DO
! CONSTRUCT THE MVL AND MVR TABLES:
      DO IVL=MVSTA,MVEND
        MVLL=IVL-MVSTA+1
        MVRR=0
        IF(IVR(IVL,1).NE.0) MVRR=IVR(IVL,1)-MVSTA+1
        MVR(MVLL,1)=MVRR
        MVRR=0
        IF(IVR(IVL,2).NE.0) MVRR=IVR(IVL,2)-MVSTA+1
        MVR(MVLL,2)=MVRR
        MVL(MVLL,1)=0
        MVL(MVLL,2)=0
      END DO
      DO MV=1,NMIDV
        IF(MVR(MV,1).NE.0) MVL(MVR(MV,1),1)=MV
        IF(MVR(MV,2).NE.0) MVL(MVR(MV,2),2)=MV
      END DO

#ifdef _DEBUGPRINT_
      WRITE(6,*)
      WRITE(6,*)' MIDVERT PAIR TABLES MVL,MVR IN MKSEG:'
      WRITE(6,*)' MVL TABLE:'
      WRITE(6,1234)(MV,MVL(MV,1),MVL(MV,2),MV=1,NMIDV)
      WRITE(6,*)' MVR TABLE:'
      WRITE(6,1234)(MV,MVR(MV,1),MVR(MV,2),MV=1,NMIDV)
 1234 FORMAT(3(3(1X,I4),4X))
      WRITE(6,*)
      WRITE(6,*)' VERTEX PAIR TABLE IVR IN MKSEG:'
      WRITE(6,1234)(IVL,IVR(IVL,1),IVR(IVL,2),IVL=1,NVERT)
#endif

! INITIALIZE SEGMENT TABLES, AND MARK VERTICES AS UNUSABLE:
      DO IVLT=1,NVERT
        DO ISGT=1,26
          ISGMNT(IVLT,ISGT)=0
          VSGMNT(IVLT,ISGT)=ZERO
        END DO
      END DO
      DO IVLT=1,NVERT
        DO ISGT=1,26
          ITT=ITVPT(ISGT)
          IVRT=IVLT
          IF((ITT.EQ.1).OR.(ITT.EQ.2)) IVRT=IVR(IVLT,ITT)
          IF(IVRT.EQ.0) cycle
          IVLB=IDOWN(IVLT,IC1(ISGT))
          IF(IVLB.EQ.0) cycle
          IVRB=IDOWN(IVRT,IC2(ISGT))
          IF(IVRB.EQ.0) cycle
! SEGMENT IS NOW ACCEPTED AS POSSIBLE.

          ISGMNT(IVLT,ISGT)=IVLB
          IB=IDRT(IVLT,IBTAB)
          Select Case (ISVC(ISGT))
             Case(1)
                V=ONE
             Case(2)
                V=-ONE
             Case(3)
                V=ONE/DBLE(1+IB)
             Case(4)
                V=-ONE/DBLE(1+IB)
             Case(5)
                V=SQRT(DBLE(IB)/DBLE(1+IB))
             Case(6)
                V=SQRT(DBLE(2+IB)/DBLE(1+IB))
             Case(7)
                V=SQRT(DBLE(IB*(2+IB)))/DBLE(1+IB)
             Case Default
                V=0.0D0 ! Dummy assigment
                Call Abend()
          End Select
          VSGMNT(IVLT,ISGT)=V
        END DO
      END DO

#ifdef _DEBUGPRINT_
        WRITE(6,*)' SEGMENT TABLE IN MKSEG.'
        WRITE(6,*)' VLT SGT ICL ICR VLB       SEGMENT TYPE         FORMULA'
        DO IV=1,NVERT
          DO ISGT=1,26
            ID=ISGMNT(IV,ISGT)
            IF(ID.EQ.0) cycle
            ICL=IC1(ISGT)
            ICR=IC2(ISGT)
            IF(ISGT.LE.4) TEXT='  WALK SECTION.'
            IF((ISGT.GE.5).AND.(ISGT.LE.8)) TEXT=' TOP SEGMENT.'
            IF((ISGT.GE.9).AND.(ISGT.LE.18)) TEXT=' MID-SEGMENT.'
            IF((ISGT.GE.19).AND.(ISGT.LE.22)) TEXT=' BOTTOM SEGMENT.'
            IF(ISGT.GT.22) TEXT=' DOWN-WALK SECTION.'
            WRITE(6,2345) IV,ISGT,ICL,ICR,ID,TEXT,FRML(ISVC(ISGT))
2345        FORMAT(1X,5I4,5X,A20,5X,A20)
          END DO
        END DO
#endif

      End Associate

      END
