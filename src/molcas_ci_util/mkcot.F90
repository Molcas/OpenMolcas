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
      SUBROUTINE MKCOT(NSYM,nLev,SGS,CIS,ISCR)
!     PURPOSE: SET UP COUNTER AND OFFSET TABLES FOR WALKS AND CSFS
!     NOTE:    TO GET GET VARIOUS COUNTER AND OFFSET TABLES
!              THE DOWN-CHAIN TABLE IS SCANNED TO PRODUCE ALL POSSIBLE
!              WALKS. POSSIBLY, THERE ARE MORE EFFICIENT WAYS, BUT
!              SINCE ONLY UPPER AND LOWER WALKS ARE REQUIRED
!              THEIR NUMBER IS VERY LIMITTED, EVEN FOR LARGE CASES.
#ifdef _DEBUGPRINT_
      use Definitions, only: LF => u6
#endif
      use Symmetry_Info, only: Mul
      use stdalloc, only: mma_allocate
      use struct, only: SGStruct, CIStruct
      IMPLICIT None
      Integer, Intent(In):: nSym, nLev
      Type (SGStruct) SGS
      Type (CIStruct) CIS
      Integer, Intent(InOut):: ISCR(3,0:NLEV)
!
      Integer, PARAMETER :: IVERT=1, ISYM=2, ISTEP=3
      Integer IHALF, IVTSTA, IVTEND, LEV1, LEV2, IVTOP, LEV, ILND, IS, ISML, ISTP, &
              ISYDWN, ISYTOT, ISYUP, IVB, IVT, IWSYM, MV, N, NUW
#ifdef _DEBUGPRINT_
      Integer NLW
#endif
      Logical Found

      CIS%NIPWLK=1+(SGS%MIDLEV-1)/15
      CIS%NIPWLK=MAX(CIS%NIPWLK,1+(SGS%NLEV-SGS%MIDLEV-1)/15)
      CALL mma_allocate(CIS%NOW,2,NSYM,CIS%NMIDV,Label='CIS%NOW')
      CALL mma_allocate(CIS%IOW,2,NSYM,CIS%NMIDV,Label='CIS%IOW')
      CALL mma_allocate(CIS%NOCSF,NSYM,CIS%NMIDV,NSYM,Label='CIS%NOCSF')
      CALL mma_allocate(CIS%IOCSF,NSYM,CIS%NMIDV,NSYM,Label='CIS%IOCSF')
      Call mma_allocate(CIS%NCSF,nSym,Label='CIS%NCSF')

      Associate(nVert=>SGS%nVert, MidLev=>SGS%MidLev, nMidV=>CIS%nMidV, MVSta=>SGS%MVSta, &
                MVEnd=>SGS%MVEnd, nWalk=>CIS%nWalk, nIpWlk=>CIS%nIpWlk, &
                ISm=>SGS%Ism, iDown=>SGS%DOwn, NOW=>CIS%NOW, IOW=>CIS%IOW, &
                NOCSF=>CIS%NOCSF, IOCSF=>CIS%IOCSF, NCSF=>CIS%NCSF)
!
!     CLEAR ARRAYS IOW AND NOW
!
      NOW(:,:,:)=0
      IOW(:,:,:)=0
!
!     CLEAR ARRAYS IOCSF AND NOCSF
!
      IOCSF(:,:,:)=0
      NOCSF(:,:,:)=0
!
!     START MAIN LOOP OVER UPPER AND LOWER WALKS, RESPECTIVELY.
!
      DO IHALF=1,2
        IF(IHALF==1) THEN
          IVTSTA=1
          IVTEND=1
          LEV1=NLEV
          LEV2=MIDLEV
        ELSE
          IVTSTA=MVSta
          IVTEND=MVEnd
          LEV1=MIDLEV
          LEV2=0
        END IF
!
!     LOOP OVER VERTICES STARTING AT TOP OF SUBGRAPH
!
        DO IVTOP=IVTSTA,IVTEND
!     SET CURRENT LEVEL=TOP LEVEL OF SUBGRAPH
          LEV=LEV1
          ISCR(IVERT,LEV)=IVTOP
          ISCR(ISYM,LEV)=1
          ISCR(ISTEP,LEV)=-1
          DO WHILE(LEV<=LEV1)
!     FIND FIRST POSSIBLE UNTRIED ARC DOWN FROM CURRENT VERTEX
          IVT=ISCR(IVERT,LEV)
          Found=.FALSE.
          DO ISTP=ISCR(ISTEP,LEV)+1,3
            IVB=IDOWN(IVT,ISTP)
            IF(IVB.NE.0) Then
              FOUND=.TRUE.
              EXIT
            END IF
          END DO
!     NO SUCH ARC WAS POSSIBLE. GO UP ONE STEP AND TRY AGAIN.
          IF(.NOT.Found) THEN
          ISCR(ISTEP,LEV)=-1
          LEV=LEV+1
          Cycle
          END IF
!     SUCH AN ARC WAS FOUND. WALK DOWN:
          ISCR(ISTEP,LEV)=ISTP
          ISML=1
          IF((ISTP==1).OR.(ISTP==2)) ISML=ISM(LEV)
          LEV=LEV-1
          ISCR(ISYM,LEV)=MUL(ISML,ISCR(ISYM,LEV+1))
          ISCR(IVERT,LEV)=IVB
          ISCR(ISTEP,LEV)=-1
          IF (LEV>LEV2) Cycle
!     WE HAVE REACHED THE BOTTOM LEVEL. THE WALK IS COMPLETE.
!     FIND MIDVERTEX NUMBER ORDERING NUMBER AND SYMMETRY OF THIS WALK
          MV=ISCR(IVERT,MIDLEV)+1-MVSta
          IWSYM=ISCR(ISYM,LEV2)
          ILND=1+NOW(IHALF,IWSYM,MV)
!     SAVE THE MAX WALK NUMBER FOR GIVEN SYMMETRY AND MIDVERTEX
          NOW(IHALF,IWSYM,MV)=ILND
!     BACK UP ONE LEVEL AND TRY AGAIN:
          LEV=LEV+1
          END DO
        END DO
      END DO
!
!     NOW,CONSTRUCT OFFSET TABLES FOR UPPER AND LOWER WALKS
!     SEPARATED FOR EACH MIDVERTEX AND SYMMETRY
!
      NUW=0
      DO MV=1,NMIDV
        DO IS=1,NSYM
          IOW(1,IS,MV)=NUW*NIPWLK
          NUW=NUW+NOW(1,IS,MV)
        END DO
      END DO
      NWALK=NUW
      DO MV=1,NMIDV
        DO IS=1,NSYM
          IOW(2,IS,MV)=NWALK*NIPWLK
          NWALK=NWALK+NOW(2,IS,MV)
        END DO
      END DO
!
!     FINALLY, CONSTRUCT COUNTER AND OFFSET TABLES FOR THE CSFS
!     SEPARATED BY MIDVERTICES AND SYMMETRY.
!     FORM ALSO CONTRACTED SUMS OVER MIDVERTICES.
!
      NCSF(:)=0
      DO ISYTOT=1,NSYM
        DO MV=1,NMIDV
          DO ISYUP=1,NSYM
            ISYDWN=MUL(ISYTOT,ISYUP)
            N=NOW(1,ISYUP,MV)*NOW(2,ISYDWN,MV)
            NOCSF(ISYUP,MV,ISYTOT)=N
            IOCSF(ISYUP,MV,ISYTOT)=NCSF(ISYTOT)
            NCSF(ISYTOT)=NCSF(ISYTOT)+N
          END DO
        END DO
      END DO
#ifdef _DEBUGPRINT_
      NLW=NWALK-NUW
      Write(LF,*)
      Write(LF,*)' TOTAL NR OF WALKS: UPPER ',NUW
      Write(LF,*)'                    LOWER ',NLW
      Write(LF,*)'                     SUM  ',NWALK
      Write(LF,*)
      Write(LF,*)' NR OF CONFIGURATIONS/SYMM:'
      Write(LF,'(8(1X,I8))')(NCSF(IS),IS=1,NSYM)
      Write(LF,*)
      Write(LF,*)
      Write(LF,*)' NR OF WALKS AND CONFIGURATIONS IN NRCOUP'
      Write(LF,*)' BY MIDVERTEX AND SYMMETRY.'
      DO MV=1,NMIDV
        Write(LF,*)
        Write(LF,'(A,I2,A,8I6)') '  MV=',MV,'    UPPER WALKS:', &
                                 (NOW(1,IS,MV),IS=1,NSYM)
        Write(LF,'(A,8I6)') '           LOWER WALKS:', &
                                 (NOW(2,IS,MV),IS=1,NSYM)
        DO ISTP=1,NSYM
        Write(LF,'(A,I2,A,8I6)') ' ISTP=',ISTP,'  CONFIGURATIONS:', &
                               (NOCSF(IS,MV,ISTP),IS=1,NSYM)
        END DO
      END DO
#endif

      End Associate

      END SUBROUTINE MKCOT
