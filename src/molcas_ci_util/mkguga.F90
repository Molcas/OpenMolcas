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
      SUBROUTINE MKGUGA(NLEV,NSYM,STSYM,Skip_MKSGNUM)
!
!     PURPOSE: MAKE THE GUGA TABLES
!     NOTE:    TO RETAIN THE TABLES AVAILABLE FOR LATER PURPOSES
!              THE START ADRESSES OF OF THE ARRAYS ETC. ARE STORED IN
!              THREE USER DEFINED TYPES. Consult the struct.F90 and gugx.F90 files for the details.
!
#ifdef _DEBUGPRINT_
      use Definitions, only: LF => u6
#endif
      use stdalloc, only: mma_allocate, mma_deallocate
      use gugx, only: IFCAS, CIS, SGS, EXS

      IMPLICIT None
!
      Integer NLEV, NSYM, STSYM
     Logical, Optional:: Skip_MKSGNUM

      Integer, Pointer:: DRTP(:,:)=>Null(), DOWNP(:)=>Null()
      Integer, Allocatable, Target:: DRT0(:,:), DOWN0(:)
      Integer, Allocatable:: TMP(:), V11(:), SCR(:)
      Integer IAC, NDOWN0, NDRT0, NLSGN, NSCR, NTMP, NUSGN, NDOWN, NDRT,  &
     &        NICASE, NVERT0

!Note that we do not associate the arrays here since the are not allocated yet.
      Associate (nVert => SGS%nVert, MidLev => SGS%MidLev, MVSta => SGS%MVSta,  &
     &           MVEnd => SGS%MVEnd,  nMidV => CIS%nMidV, nIpWlk => CIS%nIpWlk, &
     &           nWalk => CIS%nWalk,  MxUp => SGS%MxUp, MxDwn => SGS%MxDwn, &
     &            LM1RAS=>SGS%LM1RAS, LM3RAS=>SGS%LM3RAS,               &
     &            LV1RAS=>SGS%LV1RAS, LV3RAS=>SGS%LV3RAS,               &
     &            IA0 => SGS%IA0, IB0 => SGS%IB0, IC0 => SGS%IC0)
!
!     SET UP A FULL PALDUS DRT TABLE:
!     (INITIALLY NO RESTRICTIONS ARE PUT UP)
!
      IAC=MIN(IA0,IC0)
      NVERT0=((IA0+1)*(IC0+1)*(2*IB0+IAC+2))/2-(IAC*(IAC+1)*(IAC+2))/6
      NDRT0=5*NVERT0
      NDOWN0=4*NVERT0
      NTMP=((NLEV+1)*(NLEV+2))/2

      IF(IFCAS.NE.0) THEN
         CALL mma_allocate(DRT0,NVERT0,5,Label='DRT0')
         CALL mma_allocate(DOWN0,NDOWN0,Label='DOWN0')
         DRTP => DRT0
         DOWNP=> DOWN0
      ELSE
         NVERT=NVERT0
         NDRT=NDRT0
         NDOWN=NDOWN0
         CALL mma_allocate(SGS%DRT,NVERT,5,Label='SGS%DRT')
         CALL mma_allocate(SGS%DOWN,NDOWN,Label='SGS%DOWN')
         DRTP => SGS%DRT
         DOWNP=> SGS%DOWN
      ENDIF
      CALL mma_allocate(TMP,NTMP,Label='TMP')
      CALL mkDRT0 (IA0,IB0,IC0,NVERT0,DRTP,DOWNP,NTMP,TMP)
      CALL mma_deallocate(TMP)
!
#ifdef _DEBUGPRINT_
      Write(LF,*)
      Write(LF,*)' PALDUS DRT TABLE (UNRESTRICTED):'
      CALL PRDRT(NVERT0,DRTP,DOWNP)
#endif
      DOWNP=> Null()
      DRTP => Null()
!
!     IF THIS IS A RAS CALCULATION PUT UP RESTRICTIONS BY DELETING
!     VERTICES WHICH VIOLATE THE FORMER.
!
      IF(IFCAS.NE.0) THEN
        CALL mma_allocate(V11,NVERT0,Label='V11')
        CALL RESTR(NVERT0,DRT0,DOWN0,V11,LV1RAS, LV3RAS, LM1RAS, LM3RAS, NVERT)
!
!     REASSEMBLE THE DRT TABLE (REMOVE DISCONNECTED VERTICES)
!
        CALL mma_allocate(SGS%DRT,nVert,5,Label='DRT')
        CALL mma_allocate(SGS%DOWN,4*nVert,Label='SGS%DOWN')
        CALL mkDRT(NVERT0,NVERT,DRT0,DOWN0,V11,SGS%DRT,SGS%DOWN)
        CALL mma_deallocate(V11)
        CALL mma_deallocate(DRT0)
        CALL mma_deallocate(DOWN0)
!
#ifdef _DEBUGPRINT_
        Write(LF,*)
        Write(LF,*)' PALDUS DRT TABLE (RESTRICTED):'
        CALL PRDRT(NVERT,SGS%DRT,SGS%DOWN)
#endif

!     IF THIS IS A CAS CALCULATION PROCEED WITH THE UNRESTRICTED
!     DRT TABLE
!
      ENDIF
!
!     CALCULATE ARC WEIGHT.
!
      CALL mma_allocate(SGS%DAW,5*nVert,Label='SGS%DAW')
      CALL MKDAW(NVERT,SGS%DOWN,SGS%DAW)
!
!     COMPUTE UPCHAIN TABLE AND REVERSE ARC WEIGHTS
!
      CALL mma_allocate(SGS%UP,4*nVert,Label='SGS%UP')
      CALL mma_allocate(SGS%RAW,5*nVert,Label='SGS%RAW')
      CALL MKRAW(NVERT,SGS%DOWN,SGS%UP,SGS%RAW)
!
!     COMPUTE LTV TABLES.
!
      CALL mma_allocate(SGS%LTV,NLEV+2,Label='LTV')
      CALL MKLTV(NVERT,NLEV,SGS%DRT,SGS%LTV)
!
!     COMPUTE MIDLEVEL AND LIMITS ON MIDVERTICE.
!
      CALL MKMID(NVERT,NLEV,SGS%DAW,SGS%RAW,SGS%LTV,MIDLEV, NMIDV, MVSta, MVEnd, MXUP, MXDWN)

      CIS%nMidV =nMidV
!
!     FORM VARIOUS OFFSET TABLES:
!     NOTE: NIPWLK AND DOWNWLK ARE THE NUMER OF INTEGER WORDS USED
!           TO STORE THE UPPER AND LOWER WALKS IN PACKED FORM.
!
      NIPWLK=1+(MIDLEV-1)/15
      NIPWLK=MAX(NIPWLK,1+(NLEV-MIDLEV-1)/15)
      CIS%nIpWlk=nIpWlk
      NSCR=MAX(6,3*(NLEV+1))
      CALL mma_allocate(CIS%NOW,2*NMIDV*NSYM,Label='CIS%NOW')
      CALL mma_allocate(CIS%IOW,2*NMIDV*NSYM,Label='CIS%IOW')
      CALL mma_allocate(CIS%NOCSF,NMIDV*(NSYM**2),Label='CIS%NOCSF')
      CALL mma_allocate(CIS%IOCSF,NMIDV*(NSYM**2),Label='CIS%IOCSF')
      CALL mma_allocate(SCR,NSCR,Label='SCR')
      Call mma_allocate(CIS%NCSF,nSym,Label='CIS%NCSF')
      CALL MKCOT(NSYM,NLEV,NVERT,MIDLEV,NMIDV,MVSta,MVEnd,NWALK,NIPWLK,SGS%ISM,SGS%DOWN,CIS%NOW,CIS%IOW,CIS%NCSF, &
                 CIS%IOCSF,CIS%NOCSF,SCR)
!
!     CONSTRUCT THE CASE LIST
!
      NICASE=NWALK*NIPWLK
      CALL mma_allocate(CIS%ICASE,nWalk*nIpWlk,Label='CIS%ICASE')
      Call MKCLIST(NSYM,NLEV,NVERT,MIDLEV,MVSta,MVEnd,NMIDV,NICASE,NIPWLK,SGS%ISM,SGS%DOWN,CIS%NOW,CIS%IOW,CIS%ICASE,SCR)
      CALL mma_deallocate(SCR)
!
!     SET UP ENUMERATION TABLES
!
      If (Present(Skip_MKSGNUM)) Then
         If (Skip_MKSGNUM) Return
      End If
      NUSGN=MXUP*NMIDV
      NLSGN=MXDWN*NMIDV
      CALL mma_allocate(EXS%USGN,NUSGN,Label='EXS%USGN')
      CALL mma_allocate(EXS%LSGN,NLSGN,Label='EXS%LSGN')
      Call MKSGNUM(STSYM,NSYM,NLEV,NVERT,MIDLEV,NMIDV,MXUP,MXDWN,NICASE,NIPWLK,SGS%DOWN,SGS%UP,SGS%DAW,SGS%RAW,CIS%NOW, &
                   CIS%IOW,EXS%USGN,EXS%LSGN,CIS%ICASE)
!
!     EXIT
!
      End Associate

      END SUBROUTINE MKGUGA

      SUBROUTINE MKGUGA_FREE()
!
!     PURPOSE: FREE THE GUGA TABLES
!
      use stdalloc, only: mma_deallocate
      use gugx, only:  SGS, CIS, EXS
      IMPLICIT None
!
      If (Allocated(SGS%ISM)) Call mma_deallocate(SGS%ISM)
      If (Allocated(SGS%DRT)) Call mma_deallocate(SGS%DRT)
      If (Allocated(SGS%DOWN)) Call mma_deallocate(SGS%DOWN)
      If (Allocated(SGS%UP)) Call mma_deallocate(SGS%UP)
      If (Allocated(SGS%LTV)) Call mma_deallocate(SGS%LTV)

      If (Allocated(SGS%DAW)) Call mma_deallocate(SGS%DAW)
      If (Allocated(SGS%RAW)) Call mma_deallocate(SGS%RAW)

      If (Allocated(CIS%NOW)) Call mma_deallocate(CIS%NOW)
      If (Allocated(CIS%IOW)) Call mma_deallocate(CIS%IOW)
      If (Allocated(CIS%NOCSF)) Call mma_deallocate(CIS%NOCSF)
      If (Allocated(CIS%IOCSF)) Call mma_deallocate(CIS%IOCSF)
      If (Allocated(CIS%NCSF)) Call mma_deallocate(CIS%NCSF)
      If (Allocated(CIS%ICASE)) Call mma_deallocate(CIS%ICASE)

      If (Allocated(EXS%USGN)) Call mma_deallocate(EXS%USGN)
      If (Allocated(EXS%LSGN)) Call mma_deallocate(EXS%LSGN)

      If (Allocated(EXS%NOCP)) Call mma_deallocate(EXS%NOCP)
      If (Allocated(EXS%IOCP)) Call mma_deallocate(EXS%IOCP)
      If (Allocated(EXS%ICOUP)) Call mma_deallocate(EXS%ICOUP)

      If (Allocated(EXS%MVL)) Call mma_deallocate(EXS%MVL)
      If (Allocated(EXS%MVR)) Call mma_deallocate(EXS%MVR)

      If (Allocated(EXS%VTAB)) Call mma_deallocate(EXS%VTAB)
      If (Allocated(EXS%SGTMP)) Call mma_deallocate(EXS%SGTMP)

      END SUBROUTINE MKGUGA_FREE
