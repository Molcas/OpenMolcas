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
!#define _DEBUGPRINT_
      SUBROUTINE MKSGNUM(STSYM,SGS,CIS,EXS)

!     PURPOSE: FOR ALL UPPER AND LOWER WALKS
!              COMPUTE THE DIRECT ARC WEIGHT SUM AND THE
!              REVERSE ARC WEIGHT SUM, RESPECTIVELY.
!              STORE THE DATA IN THE TABLES IUSGNUM AND ILSGNUM
!
#ifdef _DEBUGPRINT_
      use Definitions, only: LF => u6
#endif
      use Symmetry_Info, only: MUL
      use struct, only: SGStruct, CIStruct, EXStruct
      use stdalloc, only: mma_allocate, mma_deallocate
      IMPLICIT None
      Integer STSYM
      Type (SGStruct) SGS
      Type (CIStruct) CIS
      Type (EXStruct) EXS
!
! local stuff
      Integer, Allocatable:: ISTEPVEC(:)
      Integer :: IC, ICODE, ICONF, IDAWSUM, ILOFF, ILW, IPOS, IRAWSUM,  &
     &           ISTEP, ISYM, IUOFF, IUW, JPOS, JSYM, LEV, LV, MIDV,    &
     &           NLW, NUW
#ifdef _DEBUGPRINT_
      Integer :: J
#endif

      CALL mma_allocate(EXS%USGN,SGS%MXUP,CIS%NMIDV,Label='EXS%USGN')
      CALL mma_allocate(EXS%LSGN,SGS%MXDWN,CIS%NMIDV,Label='EXS%LSGN')
      Call mma_allocate(ISTEPVEC,SGS%nLev,Label='ISTEPVEC')

      Associate (nSym=>SGS%nSym, nLev=>SGS%nLev, MidLev=>SGS%MidLev,    &
     &           nMidV=>CIS%nMidV, MxUp=>SGS%MxUp, MxDwn=>SGS%MxDwn,    &
     &           nIpWlk=>CIS%nIpWlk, iDOwn=>SGS%Down, iUp=>SGS%Up,      &
     &           iDaw=>SGS%Daw, iRaw=>SGS%Raw, NOW=>CIS%NOW,            &
     &           IOW=>CIS%IOW, IUSGNUM=>EXS%USGN, ILSGNUM=>EXS%LSGN,    &
     &           iCASE=>CIS%iCase,nVert=>SGS%nVert)

!
!     INITIALIZE NUMBERING TABLES
!
      DO MIDV=1,NMIDV
        DO IUW=1,MXUP
           IUSGNUM(IUW,MIDV)=0
        END DO
        DO ILW=1,MXDWN
           ILSGNUM(ILW,MIDV)=0
        END DO
      END DO
!
!     MAIN LOOP RUNS OVER MIDVERTICES AND SYMMETRIES
!
      ICONF=0
      DO MIDV=1,NMIDV
        DO ISYM=1,NSYM
          IUOFF=1+IOW(1,ISYM,MIDV)
          NUW=NOW(1,ISYM,MIDV)
          JSYM=MUL(ISYM,STSYM)
          ILOFF=1+IOW(2,JSYM,MIDV)
          NLW=NOW(2,JSYM,MIDV)
          IF( NUW.EQ.0 .OR. NLW.EQ.0 ) Cycle
!
!         LOOP OVER ALL UPPER WALKS
!
          DO IUW=1,NUW
            IPOS=IUOFF+NIPWLK*(IUW-1)
!     UNPACK THE UPPER WALK STEP VECTOR
            ICODE=ICASE(IPOS)
            JPOS=0
            DO LEV=(MIDLEV+1),NLEV
              JPOS=JPOS+1
              IF( JPOS.EQ.16 ) THEN
                JPOS=1
                IPOS=IPOS+1
                ICODE=ICASE(IPOS)
              ENDIF
              ISTEP=MOD(ICODE,4)
              ISTEPVEC(LEV)=ISTEP
              ICODE=ICODE/4
            END DO
!     GET REVERSE ARC WEIGHT FOR UPPER WALK
            IRAWSUM=1
            LV=1
            DO LEV=NLEV,(MIDLEV+1),-1
              IC=ISTEPVEC(LEV)
              LV=IDOWN(LV,IC)
              IRAWSUM=IRAWSUM+IRAW(LV,IC)
            END DO
            IUSGNUM(IRAWSUM,MIDV)=IUW
          END DO
!
!         LOOP OVER ALL LOWER WALKS
!
          DO ILW=1,NLW
            IPOS=ILOFF+NIPWLK*(ILW-1)
!     UNPACK WALK STEP VECTOR
            ICODE=ICASE(IPOS)
            JPOS=0
            DO LEV=1,MIDLEV
              JPOS=JPOS+1
              IF( JPOS.EQ.16 ) THEN
                JPOS=1
                IPOS=IPOS+1
                ICODE=ICASE(IPOS)
              ENDIF
              ISTEP=MOD(ICODE,4)
              ISTEPVEC(LEV)=ISTEP
              ICODE=ICODE/4
            END DO
!     GET DIRECT ARC WEIGHT FOR THE LOWER WALK
            IDAWSUM=1
            LV=NVERT
            DO LEV=1,MIDLEV
              IC=ISTEPVEC(LEV)
              LV=IUP(LV,IC)
              IDAWSUM=IDAWSUM+IDAW(LV,IC)
            END DO
            ILSGNUM(IDAWSUM,MIDV)=ICONF
            ICONF=ICONF+NUW
          END DO
!
        END DO
      END DO
#ifdef _DEBUGPRINT_
      Write(LF,*)
      Write(LF,*)' ILSGNUM IN SUBROUTINE MKSGNUM'
      DO MIDV=1,NMIDV
        Write(LF,'(1X,''MIDV='',I3,/,(20I6))')MIDV,(ILSGNUM(J,MIDV),J=1,MXDWN)
      END DO
      Write(LF,*)
      Write(LF,*)' IUSGNUM IN SUBROUTINE MKSGNUM'
      DO MIDV=1,NMIDV
        Write(LF,'(1X,''MIDV='',I3,/,(20I6))')MIDV,(IUSGNUM(J,MIDV),J=1,MXUP)
      END DO
      Write(LF,*)
#endif

      End Associate

      Call mma_deallocate(ISTEPVEC)

      END SUBROUTINE MKSGNUM
