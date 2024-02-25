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
      SUBROUTINE MKGUGA(NLEV,NSYM,STSYM,NCSF,Skip_MKSGNUM)
!
!     PURPOSE: MAKE THE GUGA TABLES
!     NOTE:    TO RETAIN THE TABLES AVAILABLE FOR LATER PURPOSES
!              THE START ADRESSES OF OF THE ARRAYS ARE STORED IN
!              THE COMMON /GUGA/. THESE ARE:
!              SGS%DRT,DOWN,DAW,UP,RAW,NOW1,IOW1,NOCSF,IOCSF
!
#ifdef _DEBUGPRINT_
      use Definitions, only: LF => u6
#endif
      use stdalloc, only: mma_allocate, mma_deallocate
      use gugx, only:       IA0, IB0, IC0, NVERT0,                      &
     &                IFCAS,                               &
     &                                     RAW, NRAW,               &
     &                NMIDV, MXUP, MXDWN, NWALK, NNOW,  DAW, NDAW,      &
     &                NIOW, NIPWLK, NICASE,  ICASE,       NNOCSF,       &
     &                NOCSF, NIOCSF,  IOCSF,  LSGN,  USGN, NOW1,        &
     &                IOW1, LV1RAS, LV3RAS, LM1RAS, LM3RAS,             &
     &                SGS

      IMPLICIT None
!
      Integer NLEV, NSYM, STSYM
      Integer NCSF(NSYM)
      Logical, Optional:: Skip_MKSGNUM

      Integer, Pointer:: DRTP(:)=>Null(), DOWNP(:)=>Null()
      Integer, Allocatable, Target:: DRT0(:), DOWN0(:)
      Integer, Allocatable:: TMP(:), V11(:), SCR(:)
      Integer IAC, NDOWN0, NDRT0, NLSGN,       NSCR, NTMP, NUSGN, NDOWN, NDRT
      Integer nVert, MidLev, MVSta, MVEnd
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
         CALL mma_allocate(DRT0,NDRT0,Label='DRT0')
         CALL mma_allocate(DOWN0,NDOWN0,Label='DOWN0')
         DRTP => DRT0
         DOWNP=> DOWN0
      ELSE
         NDRT=NDRT0
         NDOWN=NDOWN0
         CALL mma_allocate(SGS%DRT,NDRT,Label='SGS%DRT')
         CALL mma_allocate(SGS%DOWN,NDOWN,Label='SGS%DOWN')
         DRTP => SGS%DRT
         DOWNP=> SGS%DOWN
         NVERT=NVERT0
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
        CALL mma_allocate(SGS%DRT,5*nVert,Label='DRT')
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
      NDAW=5*NVERT
      CALL mma_allocate(DAW,NDAW,Label='DAW')
      CALL MKDAW(NVERT,SGS%DOWN,DAW)
!
!     COMPUTE UPCHAIN TABLE AND REVERSE ARC WEIGHTS
!
      NRAW=5*NVERT
      CALL mma_allocate(SGS%UP,4*nVert,Label='SGS%UP')
      CALL mma_allocate(RAW,NRAW,Label='RAW')
      CALL MKRAW(NVERT,SGS%DOWN,SGS%UP,RAW)
!
!     COMPUTE LTV TABLES.
!
      CALL mma_allocate(SGS%LTV,NLEV+2,Label='LTV')
      CALL MKLTV(NVERT,NLEV,SGS%DRT,SGS%LTV)
!
!     COMPUTE MIDLEVEL AND LIMITS ON MIDVERTICE.
!
      CALL MKMID(NVERT,NLEV,DAW,RAW,SGS%LTV,MIDLEV, NMIDV, MVSta, MVEnd, MXUP, MXDWN)

      SGS%nVert =nVert
      SGS%MidLev=MidLev
      SGS%MVSta =MVSta
      SGS%MVEnd =MVEnd
!
!     FORM VARIOUS OFFSET TABLES:
!     NOTE: NIPWLK AND DOWNWLK ARE THE NUMER OF INTEGER WORDS USED
!           TO STORE THE UPPER AND LOWER WALKS IN PACKED FORM.
!
      NIPWLK=1+(MIDLEV-1)/15
      NIPWLK=MAX(NIPWLK,1+(NLEV-MIDLEV-1)/15)
      NNOW=2*NMIDV*NSYM
      NIOW=NNOW
      NNOCSF=NMIDV*(NSYM**2)
      NIOCSF=NNOCSF
      NSCR=MAX(6,3*(NLEV+1))
      CALL mma_allocate(NOW1,NNOW,Label='NOW1')
      CALL mma_allocate(IOW1,NIOW,Label='IOW1')
      CALL mma_allocate(NOCSF,NNOCSF,Label='NOCSF')
      CALL mma_allocate(IOCSF,NIOCSF,Label='IOCSF')
      CALL mma_allocate(SCR,NSCR,Label='SCR')
      CALL MKCOT(NSYM,NLEV,NVERT,MIDLEV,NMIDV,MVSta,MVEnd,NWALK,NIPWLK,SGS%ISM,SGS%DOWN,NOW1,IOW1,NCSF,IOCSF,NOCSF,SCR)
!
!     CONSTRUCT THE CASE LIST
!
      NICASE=NWALK*NIPWLK
      CALL mma_allocate(ICASE,NICASE,Label='ICASE')
      Call MKCLIST(NSYM,NLEV,NVERT,MIDLEV,MVSta,MVEnd,NMIDV,NICASE,NIPWLK,SGS%ISM,SGS%DOWN,NOW1,IOW1,ICASE,SCR)
      CALL mma_deallocate(SCR)
!
!     SET UP ENUMERATION TABLES
!
      If (Present(Skip_MKSGNUM)) Then
         If (Skip_MKSGNUM) Return
      End If
      NUSGN=MXUP*NMIDV
      NLSGN=MXDWN*NMIDV
      CALL mma_allocate(USGN,NUSGN,Label='USGN')
      CALL mma_allocate(LSGN,NLSGN,Label='LSGN')
      Call MKSGNUM(STSYM,NSYM,NLEV,NVERT,MIDLEV,NMIDV,MXUP,MXDWN,NICASE,NIPWLK,SGS%DOWN,SGS%UP,DAW,RAW,NOW1, &
                   IOW1,USGN,LSGN,ICASE)
!
!     EXIT
!
      END SUBROUTINE MKGUGA

      SUBROUTINE MKGUGA_FREE()
!
!     PURPOSE: FREE THE GUGA TABLES
!
      use stdalloc, only: mma_deallocate
      use gugx, only:  RAW,  DAW,  NOCSF,   &
     &                 IOCSF,  ICASE, USGN, LSGN, NOW1, IOW1, &
     &                 SGS, MVL, MVR, NOCP, IOCP, ICOUP, VTAB,&
     &                 SGTMP
      IMPLICIT None
!
      If (Allocated(SGS%ISM)) Call mma_deallocate(SGS%ISM)
      If (Allocated(SGS%DRT)) Call mma_deallocate(SGS%DRT)
      If (Allocated(SGS%DOWN)) Call mma_deallocate(SGS%DOWN)

      If (Allocated(DAW)) Call mma_deallocate(DAW)
      If (Allocated(SGS%UP)) Call mma_deallocate(SGS%UP)
      If (Allocated(RAW)) Call mma_deallocate(RAW)
      If (Allocated(SGS%LTV)) Call mma_deallocate(SGS%LTV)

      If (Allocated(NOW1)) Call mma_deallocate(NOW1)
      If (Allocated(IOW1)) Call mma_deallocate(IOW1)
      If (Allocated(NOCSF)) Call mma_deallocate(NOCSF)
      If (Allocated(IOCSF)) Call mma_deallocate(IOCSF)

      If (Allocated(ICASE)) Call mma_deallocate(ICASE)

      If (Allocated(USGN)) Call mma_deallocate(USGN)
      If (Allocated(LSGN)) Call mma_deallocate(LSGN)

      If (Allocated(MVL)) Call mma_deallocate(MVL)
      If (Allocated(MVR)) Call mma_deallocate(MVR)

      If (Allocated(NOCP)) Call mma_deallocate(NOCP)
      If (Allocated(IOCP)) Call mma_deallocate(IOCP)

      If (Allocated(ICOUP)) Call mma_deallocate(ICOUP)
      If (Allocated(VTAB)) Call mma_deallocate(VTAB)
      If (Allocated(SGTMP)) Call mma_deallocate(SGTMP)

      END SUBROUTINE MKGUGA_FREE
