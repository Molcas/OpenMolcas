!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************
      SUBROUTINE SGPRWF(SGS,CIS,LSYM,PRWTHR,iSpin,CI,lCI,KeyPRSD,LUVECDET)
!
!     PURPOSE: PRINT THE WAVEFUNCTION (SPIN COUPLING AND OCCUPATIONS)
!
!     NOTE:    THIS ROUTINE USES THE SPLIT GRAPH GUGA CONVENTION, I.E.,
!              CI BLOCKS ARE MATRICES CI(I,J), WHERE THE  FIRST INDEX
!              REFERS TO THE UPPER PART OF THE WALK.
!
      use stdalloc, only: mma_allocate, mma_deallocate
      use struct, only: SGStruct, CIStruct
      use definitions, only: u6
      IMPLICIT None
      Type(SGStruct) SGS
      Type(CIStruct) CIS
!
      Integer, Intent(In) :: lCI, LSYM, iSpin, LUVECDET
      Logical, Intent(In) :: KeyPrSD
      Real*8, Intent(In) :: PRWTHR
      Real*8 CI(lCI)

      Integer ICS(50)
      Character(LEN=400) Line
      Integer, Allocatable:: Lex(:)
      Real*8 :: COEF
      Integer :: IC1, ICDPOS, ICDWN, ICONF, ICUP, ICUPOS, IDW0, IDWN,   &
     &           ISYM, ISYUP, IUP, IUW0, LEV, MV, NCI, NDWN, NNN, NUP,  &
     &           IDWNSV, iOff, ISYDWN, IMS

!
!     RECONSTRUCT THE CASE LIST
!
      If (.NOT.Allocated(CIS%ICASE)) Call MKCLIST(SGS,CIS)

      ! scratch for determinant expansion
      IF (KeyPRSD) CALL mma_allocate(LEX,SGS%NLEV,Label='LEX')

!
      Associate (nLev=>SGS%nLev, MidLev=>SGS%MidLev,                    &
     &           nIpWlk=>CIS%nIpWlk, NOCSF=>CIS%NOCSF,                  &
     &           IOCSF=>CIS%IOCSF, NOW=>CIS%NOW, NSM=>SGS%ISm,          &
     &           IOW=>CIS%IOW, nSym=>SGS%nSym, nMidV=>CIS%nMidV,        &
     &           ICASE=>CIS%ICASE)


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
      Write(u6,'(A)') Line(1:iOff+15)
      Line=' '
!
!     ENTER THE MAIN LOOP IS OVER BLOCKS OF THE ARRAY CI
!     WITH SPECIFIED MIDVERTEX MV, AND UPPERWALK SYMMETRY ISYUP.
!
      DO MV=1,NMIDV
        DO ISYUP=1,NSYM
          NCI=NOCSF(ISYUP,MV,LSYM)
          IF(NCI.EQ.0) CYCLE
          NUP=NOW(1,ISYUP,MV)
          ISYDWN=1+IEOR(ISYUP-1,LSYM-1)
          NDWN=NOW(2,ISYDWN,MV)
          ICONF=IOCSF(ISYUP,MV,LSYM)
          IUW0=1-NIPWLK+IOW(1,ISYUP,MV)
          IDW0=1-NIPWLK+IOW(2,ISYDWN,MV)
          IDWNSV=0
          DO IDWN=1,NDWN
            DO IUP=1,NUP
              ICONF=ICONF+1
              COEF=CI(ICONF)
! -- SKIP OR PRINT IT OUT?
              IF(ABS(COEF).LT.PRWTHR) CYCLE
              IF(IDWNSV.NE.IDWN) THEN
                ICDPOS=IDW0+IDWN*NIPWLK
                ICDWN=ICASE(ICDPOS)
! -- UNPACK LOWER WALK.
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
! -- UNPACK UPPER WALK:
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
! -- PRINT IT!
              Write(Line(1:),'(I8)') iConf
              iOff=10
              iSym=nSm(1)
              Do Lev=1,nLev
                 If ( nSm(Lev).ne.iSym ) iOff=iOff+1

                 Select case (ICS(Lev))
                    Case (3)
                       Write (Line(iOff+Lev:),'(A1)') '2'
                    Case (2)
                       Write (Line(iOff+Lev:),'(A1)') 'd'
                    Case (1)
                       Write (Line(iOff+Lev:),'(A1)') 'u'
                    Case (0)
                       Write (Line(iOff+Lev:),'(A1)') '0'
                    Case Default
                       Call Abend()
                 End Select

                 If ( nSm(Lev).ne.iSym ) iSym=nSm(Lev)

              End Do
              iOff=iOff+nLev+3
              Write(Line(iOff:),'(2F8.5)') COEF,COEF**2
              Write(u6,'(6X,A)') Line(1:iOff+15)
              IF (KeyPRSD) THEN
                ! use maximum spin projection value
                IMS = ISPIN-1
                WRITE(u6,*)
                CALL EXPCSF (ICS, NLEV, IMS, LEX,coef,LuVecDet)
                WRITE(u6,*)
              ENDIF
              Line=' '
            END DO
          END DO
        END DO
      END DO

      End Associate

      ! free memory for determinant expansion
      IF (KeyPRSD) CALL mma_deallocate(LEX)

      END SUBROUTINE SGPRWF
