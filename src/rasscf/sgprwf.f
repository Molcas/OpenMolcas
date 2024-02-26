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
      SUBROUTINE SGPRWF(iSel,NOCSF,IOCSF,NOW,IOW,CI,nMidV)
C
C     PURPOSE: PRINT THE WAVEFUNCTION (SPIN COUPLING AND OCCUPATIONS)
C
C     NOTE:    THIS ROUTINE USES THE SPLIT GRAPH GUGA CONVENTION, I.E.,
C              CI BLOCKS ARE MATRICES CI(I,J), WHERE THE  FIRST INDEX
C              REFERS TO THE UPPER PART OF THE WALK.
C
      use stdalloc, only: mma_allocate, mma_deallocate
      use gugx, only: NWALK, CIS,  ICASE, SGS
      IMPLICIT REAL*8 (A-H,O-Z)
C
#include "rasdim.fh"
#include "rasscf.fh"
#include "general_mul.fh"
#include "input_ras.fh"
#include "output_ras.fh"
#include "WrkSpc.fh"
C
      Integer, Intent(In) :: nMidV
      Integer   iSel(*)
      DIMENSION NOCSF(NSYM,NMIDV,NSYM),IOCSF(NSYM,NMIDV,NSYM)
      DIMENSION NOW(2,NSYM,NMIDV),IOW(2,NSYM,NMIDV)
      DIMENSION CI(NCONF)

      DIMENSION ICS(mxact)
      Character(LEN=400) Line
      Integer nUp, nVert, MidLev, MVSta, MVEnd, nLev, nIpWlk, NICASE
      Integer, Allocatable:: Scr(:), Lex(:)
C
      NICASE=SIZE(ICASE)

      nLev  = SGS%nLev
      nVert = SGS%nVert
      MidLev= SGS%MidLev
      MVSta = SGS%MVSta
      MVEnd = SGS%MVEnd
      nIpWlk= CIS%nIpWlk

      Line(1:16)='      conf/sym  '
      iOff=16
      iSym=nSm(1)
      Do Lev=1,nLev
         If ( nSm(Lev).ne.iSym ) iOff=iOff+1
         Write(Line(iOff+Lev:),'(I1)') nSm(Lev)
         If ( nSm(Lev).ne.iSym ) iSym=nSm(Lev)
      End Do
      iOff=iOff+nLev+3
      Line(iOff:iOff+15)='   Coeff  Weight'
      Write(LF,'(A)') Line(1:iOff+15)
      Line=' '
C
C     RECONSTRUCT THE CASE LIST
C
      NSCR=3*(NLEV+1)
      NICASE=NWALK*NIPWLK
      CALL mma_allocate(SCR,NSCR,Label='SCR')
      Call MKCLIST(NSYM,NLEV,NVERT,MIDLEV,MVSta,MVEnd,NMIDV,
     &             NICASE,NIPWLK,NSM,SGS%DOWN,CIS%NOW,CIS%IOW,ICASE,SCR)

      CALL mma_deallocate(SCR)

      ! scratch for determinant expansion
      IF (KeyPRSD) CALL mma_allocate(LEX,NLEV,Label='LEX')
C
C     ENTER THE MAIN LOOP IS OVER BLOCKS OF THE ARRAY CI
C     WITH SPECIFIED MIDVERTEX MV, AND UPPERWALK SYMMETRY ISYUP.
C
      DO MV=1,NMIDV
        DO ISYUP=1,NSYM
          NCI=NOCSF(ISYUP,MV,STSYM)
          IF(NCI.EQ.0) CYCLE
          NUP=NOW(1,ISYUP,MV)
          ISYDWN=MUL(ISYUP,STSYM)
          NDWN=NOW(2,ISYDWN,MV)
          ICONF=IOCSF(ISYUP,MV,STSYM)
          IUW0=1-NIPWLK+IOW(1,ISYUP,MV)
          IDW0=1-NIPWLK+IOW(2,ISYDWN,MV)
          IDWNSV=0
          DO IDWN=1,NDWN
            DO IUP=1,NUP
              ICONF=ICONF+1
              COEF=CI(ICONF)
C -- SKIP OR PRINT IT OUT?
              IF(ABS(COEF).LT.PRWTHR) CYCLE
              IF(IDWNSV.NE.IDWN) THEN
                ICDPOS=IDW0+IDWN*NIPWLK
                ICDWN=ICASE(ICDPOS)
C -- UNPACK LOWER WALK.
                NNN=0
                DO LEV=1,MIDLEV
                  NNN=NNN+1
                  IF(NNN.EQ.16) THEN
                    NNN=1
                    ICDPOS=ICDPOS+1
                    ICDWN=ICASE(ICDPOS)
                  END IF
                  IC1=ICDWN/4
                  ICS(LEV)=ICDWN-4*IC1
                  ICDWN=IC1
                END DO
                IDWNSV=IDWN
              END IF
              ICUPOS=IUW0+NIPWLK*IUP
              ICUP=ICASE(ICUPOS)
C -- UNPACK UPPER WALK:
              NNN=0
              DO LEV=MIDLEV+1,NLEV
                NNN=NNN+1
                IF(NNN.EQ.16) THEN
                  NNN=1
                  ICUPOS=ICUPOS+1
                  ICUP=ICASE(ICUPOS)
                END IF
                IC1=ICUP/4
                ICS(LEV)=ICUP-4*IC1
                ICUP=IC1
              END DO
C -- PRINT IT!
              iSel(iConf) = 1
              Write(Line(1:),'(I8)') iConf
              iOff=10
              iSym=nSm(1)
              Do Lev=1,nLev
                 If ( nSm(Lev).ne.iSym ) iOff=iOff+1
                 If ( ICS(Lev).eq.3 ) then
                    Write(Line(iOff+Lev:),'(A1)') '2'
                 Else If (ICS(Lev).eq.2) then
                    Write(Line(iOff+Lev:),'(A1)') 'd'
                 Else If (ICS(Lev).eq.1) then
                    Write(Line(iOff+Lev:),'(A1)') 'u'
                 Else If (ICS(Lev).eq.0) then
                    Write(Line(iOff+Lev:),'(A1)') '0'
                 End If
                 If ( nSm(Lev).ne.iSym ) iSym=nSm(Lev)
              End Do
              iOff=iOff+nLev+3
              Write(Line(iOff:),'(2F8.5)') COEF,COEF**2
              Write(LF,'(6X,A)') Line(1:iOff+15)
              IF (KeyPRSD) THEN
                ! use maximum spin projection value
                IMS = ISPIN-1
                WRITE(6,*)
                CALL EXPCSF (ICS, NLEV, IMS, LEX,coef,LuVecDet)
                WRITE(6,*)
              ENDIF
              Line=' '
            END DO
          END DO
        END DO
      END DO

      ! free memory for determinant expansion
      IF (KeyPRSD) CALL mma_deallocate(LEX)

      END SUBROUTINE SGPRWF
