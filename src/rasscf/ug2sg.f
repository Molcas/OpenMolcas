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
* Copyright (C) 1990, Jeppe Olsen                                      *
*               1990, Markus P. Fuelscher                              *
************************************************************************
      SUBROUTINE UG2SG(NROOTS,NCONF,NORB,NEL,IREFSM,IPRINT,
     *                 ICONF,ISPIN,IORD,ICI,JCJ,CCI,MXROOTS)
C
C     AUTHOR:  J. OLSEN AND M.P. FUELSCHER
C              UNIV. OF LUND, SWEDEN 1990
C
C     PURPOSE: CONSTRUCT THE REINDEXING ARRAY WHICH REORDERS
C              THE CSFS GENERATED BY THE DETERMINANT CODE INTO
C              THE SPLIT GRAPH GUGA ORDER. THEREFORE, CONSTRUCT
C              FOR EACH CSF THE CORRESPONDING
C              STEP VECTOR AND PASS IT TO THE FUNCTIONS IPHASE
C              AND ISGNUM WHICH COMPUTES THE THE PHASE FACTOR
C              INVOLVED WHEN GOING FROM THE SYMMETRIC TO THE
C              UNITARY GROUP AND THE SPLIT ORDERING NUMBER.
C
      use gugx, only: SGS,CIS, EXS
      use output_ras, only: LF
      use spinfo, only: NTYP,MINOP,NCNFTP,NCSFTP
      IMPLICIT None
      INTEGER NROOTS,NCONF,NORB,NEL,IREFSM,IPRINT,MXROOTS
      INTEGER ICONF(*),ISPIN(*)
      INTEGER IORD(*)
      INTEGER ICI(MXROOTS,*),JCJ(MXROOTS,*)
      REAL*8 CCI(MXROOTS,*)
C
#include "rasdim.fh"
C
      INTEGER IWALK(mxAct)
      INTEGER KCNF(mxAct)
      Integer, External:: IPHASE
      Integer nVert, nLev, nMidV, MxUp, MxDwn
      Integer K,L,KREF,KROOT,ICSFJP,ICNBS0,IPBAS,ITYP,IOPEN,ICL,IC,
     &        ICNBS,IICSF,ICSBAS,IIBOP,IIBCL,JOCC,KOCC,ISG,IP,LPRINT,
     &        I,ISGNUM,KORB
      REAL*8 PHASE

      nLev  = SGS%nLev
      nVert = SGS%nVert
      nMidV = CIS%nMidV
      MxUp  = SGS%MxUp
      MxDwn = SGS%MxDwn
C
C     JCJ IS A TEMPORARY COPY OF ICI AND WILL OBTAIN THE SELECTED REFERENCE
C     NUMBERS IN THE SYMMETRIC GROUP NUMBERING
C

      IF( IPRINT.GE.5 ) THEN
        Write(LF,*)
        Write(LF,*)' SPLIT GRAPH GUGA CONFIGURATION NUMBERS:'
        DO 88 K=1,NROOTS
        Write(LF,'(A,I2,A,5I8)') ' ROOT',K,' CSFs:',(ICI(K,L),L=1,5)
88      CONTINUE
      ENDIF
C
      Do kRef = 1,mxRef
        Do kRoot = 1,mxRoots
          JCJ(kRoot,kRef) = 0
        End Do
      End Do
C
C     LOOP OVER CONFIGURATIONS TYPES
C
      ICSFJP = 0
      ICNBS0 = 0 ! dummy initialize
      IPBAS  = 0 ! dummy initialize
      DO 1000 ITYP = 1, NTYP
        IOPEN = ITYP + MINOP - 1
        ICL = (NEL - IOPEN) / 2
C      BASE ADDRESS FOR CONFIGURATION OF THIS TYPE
        IF( ITYP .EQ. 1 ) THEN
          ICNBS0 = 1
        ELSE
          ICNBS0 = ICNBS0 + NCNFTP(ITYP-1,IREFSM)*(NEL+IOPEN-1)/2
        END IF
C      BASE ADDRESS FOR PROTOTYPE SPIN COUPLINGS
        IF( ITYP .EQ. 1 ) THEN
          IPBAS = 1
        ELSE
          IPBAS = IPBAS + NCSFTP(ITYP-1)*(IOPEN-1)
        END IF
C
C     LOOP OVER NUMBER OF CONFIGURATIONS OF TYPE ITYP AND PROTOTYPE
C     SPIN COUPLINGS
C
        DO 900  IC = 1, NCNFTP(ITYP,IREFSM)
          ICNBS = ICNBS0 + (IC-1)*(IOPEN+ICL)
          DO 800 IICSF = 1,NCSFTP(ITYP)
            ICSFJP = ICSFJP + 1
            ICSBAS = IPBAS + (IICSF-1)*IOPEN
CPAM04 The following lines have been replaced....
*C     COMPUTE STEP VECTOR
*            CALL STEPVEC(ICONF(ICNBS),ICONF(ICNBS+ICL),ICL,IOPEN,
*     &                   ISPIN(ICSBAS),NORB,IWALK)
CPAM04 ... because of new convention for representing SG config, we
* need to arrange orbital indices in form used by RASSCF codes:
*. Obtain configuration in standard RASSCF form
            IIBOP = 1
            IIBCL = 1
            JOCC  = ICL + IOPEN
            DO KOCC = 0, JOCC-1
              KORB = ICONF(ICNBS+KOCC)
              IF(KORB.LT.0) THEN
*. Doubly occupied orbitals
                KCNF(IIBCL) = ABS(KORB)
                IIBCL = IIBCL + 1
              ELSE
*. Singly occupied orbital
                KCNF(ICL+IIBOP) = KORB
                IIBOP = IIBOP + 1
              END IF
            END DO
C     COMPUTE STEP VECTOR
            CALL STEPVEC(KCNF(1),KCNF(1+ICL),ICL,IOPEN,
     &                   ISPIN(ICSBAS),NORB,IWALK)
CPAM04 End of replacement code.
C     GET SPLIT GRAPH ORDERING NUMBER
            ISG=ISGNUM(NLEV,NVERT,SGS%MIDLEV,SGS%MVSta,NMIDV,MXUP,MXDWN,
     &                 SGS%DOWN,SGS%UP,SGS%DAW,SGS%RAW,
     &                 EXS%USGN,EXS%LSGN,IWALK)

C     GET PHASE PHASE FACTOR
            IP=IPHASE(NLEV,NVERT,SGS%DRT,SGS%UP,IWALK)
C     UPDATE REINDEXING TABLE
           IORD(ICSFJP) = ISG * IP
800       CONTINUE
900     CONTINUE
1000  CONTINUE
C
      IF( IPRINT.GE.5 ) THEN
        LPRINT=MIN(200,NCONF)
        Write(LF,*)
        Write(LF,*)' INDEX TABLE IN SUBROUTINE REORD'
        Write(LF,'(10I8)') (IORD(I),I=1,LPRINT)
        Write(LF,*)
      ENDIF
C
C     REPLACE CONFIGURATION NUMBERS
C
      DO 200 IC=1,NCONF
        ISG=IORD(IC)
        PHASE=1.0D0
        IF ( ISG.LT.0 ) PHASE=-1.0D0
        ISG=ABS(ISG)
        DO 210 K=1,NROOTS
          DO 220 L=1,MXREF
            IF (ICI(K,L).EQ.ISG) THEN
               JCJ(K,L)=IC
               CCI(K,L)=CCI(K,L)*PHASE
            END IF
220       CONTINUE
210     CONTINUE
200   CONTINUE
C
      IF( IPRINT.GE.5 ) THEN
        Write(LF,*)' SYMMETRIC GROUP CONFIGURATION NUMBERS:'
        DO 260 K=1,NROOTS
        Write(LF,'(A,I2,A,5I6)') ' ROOT',K,' CSFs:',(JCJ(K,L),L=1,5)
260     CONTINUE
        Write(LF,*)
      ENDIF
C
      RETURN
      END
