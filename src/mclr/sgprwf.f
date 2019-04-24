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
      SUBROUTINE SGPRWF_MCLR
     &           (LSYM,PRWTHR,
     &            NSYM,NLEV,NCONF,MIDLEV,NMIDV,NIPWLK,NICASE,
     &            NSM,NOCSF,IOCSF,NOW,IOW,ICASE,CI)
C
C     PURPOSE: PRINT THE WAVEFUNCTION (SPIN COUPLING AND OCCUPATIONS)
C
C     NOTE:    THIS ROUTINE USES THE SPLIT GRAPH GUGA CONVENTION, I.E.,
C              CI BLOCKS ARE MATRICES CI(I,J), WHERE THE  FIRST INDEX
C              REFERS TO THE UPPER PART OF THE WALK.
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      DIMENSION NOCSF(NSYM,NMIDV,NSYM),IOCSF(NSYM,NMIDV,NSYM)
      DIMENSION NOW(2,NSYM,NMIDV),IOW(2,NSYM,NMIDV)
      DIMENSION CI(NCONF)
      DIMENSION ICS(50)
      DIMENSION ICASE(NICASE)
      DIMENSION NSM(*)
      Character*120 Line
C
      Line(1:16)='      conf/sym  '
      iOff=16
      iSym=nSm(1)
      Do Lev=1,nLev
         If ( nSm(Lev).ne.iSym ) iOff=iOff+1
         Write (Line(iOff+Lev:),'(I1)') nSm(Lev)
         If ( nSm(Lev).ne.iSym ) iSym=nSm(Lev)
      End Do
      iOff=iOff+nLev+3
      Line(iOff:iOff+15)='   Coeff  Weight'
      Write (6,'(A)') Line(1:iOff+15)
      Line=' '
C
C
C     THE MAIN LOOP IS OVER BLOCKS OF THE ARRAY CI
C     WITH SPECIFIED MIDVERTEX MV, AND UPPERWALK SYMMETRY ISYUP.
C

      DO 40 MV=1,NMIDV
        DO 40 ISYUP=1,NSYM
          NCI=NOCSF(ISYUP,MV,LSYM)
          IF(NCI.EQ.0) GOTO 40
          NUP=NOW(1,ISYUP,MV)
          ISYDWN=1+IEOR(ISYUP-1,LSYM-1)
          NDWN=NOW(2,ISYDWN,MV)
          ICONF=IOCSF(ISYUP,MV,LSYM)
          IUW0=1-NIPWLK+IOW(1,ISYUP,MV)
          IDW0=1-NIPWLK+IOW(2,ISYDWN,MV)
          IDWNSV=0
          DO 30 IDWN=1,NDWN
            DO 30 IUP=1,NUP
              ICONF=ICONF+1
              COEF=CI(ICONF)
C -- SKIP OR PRINT IT OUT?
              IF(ABS(COEF).LT.PRWTHR) GOTO  30
              IF(IDWNSV.NE.IDWN) THEN
                ICDPOS=IDW0+IDWN*NIPWLK
                ICDWN=ICASE(ICDPOS)
C -- UNPACK LOWER WALK.
                NNN=0
                DO 10 LEV=1,MIDLEV
                  NNN=NNN+1
                  IF(NNN.EQ.16) THEN
                    NNN=1
                    ICDPOS=ICDPOS+1
                    ICDWN=ICASE(ICDPOS)
                  END IF
                  IC1=ICDWN/4
                  ICS(LEV)=ICDWN-4*IC1
                  ICDWN=IC1
10              CONTINUE
                IDWNSV=IDWN
              END IF
              ICUPOS=IUW0+NIPWLK*IUP
              ICUP=ICASE(ICUPOS)
C -- UNPACK UPPER WALK:
              NNN=0
              DO 20 LEV=MIDLEV+1,NLEV
                NNN=NNN+1
                IF(NNN.EQ.16) THEN
                  NNN=1
                  ICUPOS=ICUPOS+1
                  ICUP=ICASE(ICUPOS)
                END IF
                IC1=ICUP/4
                ICS(LEV)=ICUP-4*IC1
                ICUP=IC1
20            CONTINUE
C -- PRINT IT!
              Write (Line(1:),'(I8)') iConf
              iOff=10
              iSym=nSm(1)
              Do Lev=1,nLev
                 If ( nSm(Lev).ne.iSym ) iOff=iOff+1
                 If ( ICS(Lev).eq.3 ) then
                    Write (Line(iOff+Lev:),'(A1)') '2'
                 Else If (ICS(Lev).eq.2) then
                    Write (Line(iOff+Lev:),'(A1)') 'd'
                 Else If (ICS(Lev).eq.1) then
                    Write (Line(iOff+Lev:),'(A1)') 'u'
                 Else If (ICS(Lev).eq.0) then
                    Write (Line(iOff+Lev:),'(A1)') '0'
                 End If
                 If ( nSm(Lev).ne.iSym ) iSym=nSm(Lev)
              End Do
              iOff=iOff+nLev+3
              Write (Line(iOff:),'(2F8.5)') COEF,COEF**2
              Write (6,'(6X,A)') Line(1:iOff+15)
              Line=' '
30        CONTINUE
40    CONTINUE
C
C
C     EXIT
C
      RETURN
      END
