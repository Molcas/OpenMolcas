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
      SUBROUTINE PRWF1(ISGSTRUCT,ICISTRUCT,NLEV,NMIDV,ISM,ICS,
     &                 NOCSF,IOCSF,NOW,IOW,ISYCI,CI,CITHR)
      IMPLICIT REAL*8 (A-H,O-Z)
#include "Struct.fh"
#include "symmul.fh"
      Dimension iSGStruct(nSGSize)
      Dimension iCIStruct(nCISize)
#include "WrkSpc.fh"
      DIMENSION NOCSF(NSYM,NMIDV,NSYM),IOCSF(NSYM,NMIDV,NSYM)
      DIMENSION NOW(2,NSYM,NMIDV),IOW(2,NSYM,NMIDV)
      DIMENSION CI(*)
      DIMENSION ISM(NLEV), ICS(NLEV)
      LOGICAL SGINFO
      PARAMETER (SGINFO=.TRUE.)
      CHARACTER*80 LINE
      CHARACTER*1 CODE(0:3)
      DATA CODE /'0','u','d','2'/

C -- NOTE: THIS PRWF ROUTINE USES THE CONVENTION THAT CI BLOCKS
C -- ARE MATRICES CI(I,J), WHERE THE   F I R S T   INDEX I REFERS TO
C -- THE   U P P E R   PART OF THE WALK.
C -- THE MAIN LOOP IS OVER BLOCKS OF THE ARRAY CI
C    WITH SPECIFIED MIDVERTEX MV, AND UPPERWALK SYMMETRY ISYUP.

      MIDLEV=ISGSTRUCT(8)
      NIPWLK=ICISTRUCT(2)
      LICASE=ICISTRUCT(9)

C Size of occup/spin coupling part of line:
      WRITE(6,*)' Occupation of active orbitals, and spin coupling'
      WRITE(6,*)' of open shells. (u,d: Spin up or down).'
      WRITE(LINE,'(20A4)')('    ',I=1,20)
      K=0
      ISY=0
      DO  LEV=1,NLEV
        IF(ISY.NE.ISM(LEV)) THEN
          ISY=ISM(LEV)
          K=K+1
        END IF
        K=K+1
      END DO
      KOCLAB=10
      KOCSZ=MAX(K,KOCLAB)
      KPAD1=(KOCSZ-KOCLAB)/2
      KPAD2=(KOCSZ-K)/2
      IF(SGINFO) THEN
      WRITE(6,*)' SGUGA info is (Midvert:IsyUp:UpperWalk/LowerWalk)'
      END IF
      LINE(1:7)='  Conf '
      K=7
      IF(SGINFO) THEN
        LINE(8:22)='  SGUGA info   '
        K=22
      END IF
      LINE(K+KPAD1:K+KPAD1+9)='Occupation'
      K=K+KOCSZ
      LINE(K:K+23)='       Coef       Weight'
      WRITE(6,*) LINE
      WRITE(LINE,'(20A4)')('    ',I=1,20)

C -- THE MAIN LOOP IS OVER BLOCKS OF THE ARRAY CI
C    WITH SPECIFIED MIDVERTEX MV, AND UPPERWALK SYMMETRY ISYUP.
      DO 40 MV=1,NMIDV
        DO 40 ISYUP=1,NSYM
          NCI=NOCSF(ISYUP,MV,ISYCI)
          IF(NCI.EQ.0) GOTO 40
          NUP=NOW(1,ISYUP,MV)
          ISYDWN=MUL(ISYUP,ISYCI)
          NDWN=NOW(2,ISYDWN,MV)
          ICONF=IOCSF(ISYUP,MV,ISYCI)
          IUW0=LICASE-NIPWLK+IOW(1,ISYUP,MV)
          IDW0=LICASE-NIPWLK+IOW(2,ISYDWN,MV)
          IDWNSV=0
          DO 30 IDWN=1,NDWN
            DO 30 IUP=1,NUP
              ICONF=ICONF+1
              COEF=CI(ICONF)
C -- SKIP OR PRINT IT OUT?
              IF(ABS(COEF).LT.CITHR) GOTO  30
              IF(IDWNSV.NE.IDWN) THEN
                ICDPOS=IDW0+IDWN*NIPWLK
                ICDWN=IWORK(ICDPOS)
C -- UNPACK LOWER WALK.
                NNN=0
                DO 10 LEV=1,MIDLEV
                  NNN=NNN+1
                  IF(NNN.EQ.16) THEN
                    NNN=1
                    ICDPOS=ICDPOS+1
                    ICDWN=IWORK(ICDPOS)
                  END IF
                  IC1=ICDWN/4
                  ICS(LEV)=ICDWN-4*IC1
                  ICDWN=IC1
10              CONTINUE
                IDWNSV=IDWN
              END IF
              ICUPOS=IUW0+NIPWLK*IUP
              ICUP=IWORK(ICUPOS)
C -- UNPACK UPPER WALK:
              NNN=0
              DO 20 LEV=MIDLEV+1,NLEV
                NNN=NNN+1
                IF(NNN.EQ.16) THEN
                  NNN=1
                  ICUPOS=ICUPOS+1
                  ICUP=IWORK(ICUPOS)
                END IF
                IC1=ICUP/4
                ICS(LEV)=ICUP-4*IC1
                ICUP=IC1
20            CONTINUE
C -- PRINT IT!
              WRITE(LINE(1:7),'(I6,1X)') ICONF
              K=7
              IF(SGINFO) THEN
                LINE(K+1:K+1)='('
                WRITE(LINE(K+2:K+3),'(I2)') MV
                LINE(K+4:K+4)=':'
                WRITE(LINE(K+5:K+5),'(I1)') ISYUP
                LINE(K+6:K+6)=':'
                WRITE(LINE(K+7:K+9),'(I3)') IUP
                LINE(K+10:K+10)='/'
                WRITE(LINE(K+11:K+13),'(I3)') IDWN
                LINE(K+14:K+14)=')'
                K=K+14
              END IF
              KNXT=K+KOCSZ
              K=K+KPAD2
              ISY=0
              DO LEV=1,NLEV
                IF(ISY.NE.ISM(LEV)) THEN
                  ISY=ISM(LEV)
                  K=K+1
                  LINE(K:K)=' '
                END IF
                K=K+1
                LINE(K:K)=CODE(ICS(LEV))
              END DO
              K=KNXT
              K=K+1
              LINE(K:K+4)='     '
              K=K+5
              WRITE(LINE(K:K+7),'(F8.5)') COEF
              K=K+8
              LINE(K:K+4)='     '
              K=K+5
              WRITE(LINE(K:K+7),'(F8.5)') COEF**2
              WRITE(6,*)LINE(1:K+7)
30        CONTINUE
40    CONTINUE
      WRITE(6,*)
      WRITE(6,*)('*',I=1,80)
      RETURN
      END
