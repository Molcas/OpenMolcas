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
! Copyright (C) Per Ake Malmqvist                                      *
!               Markus P. Fuelscher                                    *
!***********************************************************************
!#define _DEBUGPRINT_
      SUBROUTINE GUGACTL()
!
!     PURPOSE: CONTROL ROUTINE TO SET UP GUGA TABLES
!     AUTHOR:  P.-AA. MALMQVIST
!
!     MODIFIED TO FIT THE DETRAS PROGRAM BY M.P. FUELSCHER
!
      use Definitions, only: LF => u6
      use gugx, only: NLEV, IA0, IB0, IC0, NVERT0, IFCAS, LV1RAS,       &
     &                LM1RAS, LV3RAS, LM3RAS, NCSF

      IMPLICIT REAL*8 (A-H,O-Z)
!
#include "rasdim.fh"
#include "warnings.h"
#include "rasscf.fh"
#include "general.fh"
#include "gas.fh"
      Character(LEN=16), Parameter :: ROUTINE='GUGACTL '

      Interface
      SUBROUTINE MKGUGA(NSM,NLEV,NSYM,STSYM,NCSF)
      IMPLICIT None

      Integer NLEV, NSYM, STSYM
      Integer NSM(NLEV)
      Integer NCSF(NSYM)
      End SUBROUTINE MKGUGA
      End Interface

! Local print level (if any)
#ifdef _DEBUGPRINT_
      WRITE(LF,*)' Entering ',ROUTINE
#endif
!
!     SET IFCAS FLAG
!     IFCAS = 0 : THIS IS A CAS CALCULATION
!     IFCAS = 1 : THIS IS A RAS CALCULATION
!
      IFCAS=0
      IF (NHOLE1.NE.0.OR.NELEC3.NE.0) IFCAS=1
      DO IS=1,NSYM
        IF (IFCAS.NE.0.AND.NASH(IS).NE.0)IFCAS=IFCAS+1
      END DO
!
!     CREATE THE SYMMETRY INDEX VECTOR
!
      CALL MKNSM()
!
!     (IFCAS-1) IS THE NUMBER OF SYMMETRIES CONTAINING ACTIVE ORBITALS
!     IF THIS IS GREATER THAN 1 ORBITAL REORDERING INTEGRALS IS REQUIRED
!     SET UP THE REINDEXING TABLE
!
      CALL SETSXCI()
!
!     FIND TOTAL NUMBER OF VERTICES IN THE SUBSPACES
!
!... for RAS
      NLEV=0
      DO IS=1,NSYM
        NLEV=NLEV+NRS1(IS)
      END DO
      LV1RAS=NLEV
      DO IS=1,NSYM
        NLEV=NLEV+NRS2(IS)
      END DO
      NRAS2=NLEV-LV1RAS
      DO IS=1,NSYM
        NLEV=NLEV+NRS3(IS)
      END DO
!
!     COMPUTE RAS RESTRICTIONS ON VERTICES:
!
      LV3RAS=LV1RAS+NRAS2
      LM1RAS=2*LV1RAS-NHOLE1
      LM3RAS=NACTEL-NELEC3
!
!     COMPUTE TOP ROW OF THE GUGA TABLE
!
      IB0=ISPIN-1
      IA0=(NACTEL-IB0)/2
      IC0=NLEV-IA0-IB0
      IF ( ((2*IA0+IB0).NE.NACTEL) .OR.                                 &
     &     (IA0.LT.0) .OR.                                              &
     &     (IB0.LT.0) .OR.                                              &
     &     (IC0.LT.0) ) then
        Write(LF,*)'GUGACTL Error: Impossible specifications.'
        Write(LF,'(1x,a,3I8)')'NACTEL,NLEV,ISPIN:',NACTEL,NLEV,ISPIN
        Write(LF,'(1x,a,3I8)')'IA0,IB0,IC0:      ',IA0,IB0,IC0
        Write(LF,*)' This is a severe internal error, or possibly'
        Write(LF,*)' indicates a strange input which should have been'
        Write(LF,*)' diagnosed earlier. Please submit a bug report.'
        Call Quit(_RC_GENERAL_ERROR_)
      End If
      IAC=MIN(IA0,IC0)
      NVERT0=((IA0+1)*(IC0+1)*(2*IB0+IAC+2))/2-(IAC*(IAC+1)*(IAC+2))/6
      If ( NVERT0.eq.0 ) then
        NCONF=0
        Return
      End If
      If ( doBlockDMRG ) then
        NCONF=1
        Return
      End If
!
!     INITIALIZE GUGA TABLES:
!
      CALL MKGUGA(NSM,NLEV,NSYM,STSYM,NCSF)
      NCONF=NCSF(STSYM)
      If ( NAC.eq.0 ) NCONF=1

#ifdef _TEST_
      Contains
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
      SUBROUTINE MKGUGA(NSM,NLEV,NSYM,STSYM,NCSF)
!
!     PURPOSE: MAKE THE GUGA TABLES
!     NOTE:    TO RETAIN THE TABLES AVAILABLE FOR LATER PURPOSES
!              THE START ADRESSES OF OF THE ARRAYS ARE STORED IN
!              THE COMMON /GUGA/. THESE ARE:
!              DRT,DOWN,DAW,UP,RAW,NOW1,IOW1,NOCSF,IOCSF
!
#ifdef _DEBUGPRINT_
      use Definitions, only: LF => u6
#endif
      use stdalloc, only: mma_allocate, mma_deallocate
      use gugx, only:       IA0, IB0, IC0, NVERT0,                      &
     &                IFCAS, NVERT, NDRT,  DRT,                         &
     &                NDOWN,  DOWN,  UP, NUP,  RAW, NRAW, MIDLEV,       &
     &                NMIDV, MXUP, MXDWN, NWALK, NNOW,  DAW, NDAW,      &
     &                NIOW, NIPWLK, NICASE,  ICASE,       NNOCSF,       &
     &                NOCSF, NIOCSF,  IOCSF,  LSGN,  USGN, NOW1,        &
     &                IOW1, LV1RAS, LV3RAS, LM1RAS, LM3RAS,             &
     &                MIDV1, MIDV2, LTV

      IMPLICIT None
!
      Integer NLEV, NSYM, STSYM
      Integer NSM(NLEV)
      Integer NCSF(NSYM)

      Integer, Pointer:: DRTP(:)=>Null(), DOWNP(:)=>Null()
      Integer, Allocatable, Target:: DRT0(:), DOWN0(:)
      Integer, Allocatable:: TMP(:), V11(:), SCR(:)
      Integer IAC, NDOWN0, NDRT0, NLSGN, NLTV, NSCR, NTMP, NUSGN
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
         CALL mma_allocate(DRT,NDRT,Label='DRT')
         CALL mma_allocate(DOWN,NDOWN,Label='DOWN')
         DRTP => DRT
         DOWNP=> DOWN
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
        NDRT=5*NVERT
        NDOWN=4*NVERT
        CALL mma_allocate(DRT,NDRT,Label='DRT')
        CALL mma_allocate(DOWN,NDOWN,Label='DOWN')
        CALL mkDRT(NVERT0,NVERT,DRT0,DOWN0,V11,DRT,DOWN)
        CALL mma_deallocate(V11)
        CALL mma_deallocate(DRT0)
        CALL mma_deallocate(DOWN0)
!
#ifdef _DEBUGPRINT_
        Write(LF,*)
        Write(LF,*)' PALDUS DRT TABLE (RESTRICTED):'
        CALL PRDRT(NVERT,DRT,DOWN)
#endif

!     IF THIS IS A CAS CALCULATION PROCEED WITH THE UNRESTRICTED
!     DRT TABLE
!
      ENDIF
!
!     CALCULATE ARC WEIGHT AND LTV TABLES.
!
      NDAW=5*NVERT
      CALL mma_allocate(DAW,NDAW,Label='DAW')
      CALL MKDAW(NVERT,DOWN,DAW)
!
!     COMPUTE UPCHAIN TABLE AND REVERSE ARC WEIGHTS
!
      NUP=4*NVERT
      NRAW=5*NVERT
      CALL mma_allocate(UP,NUP,Label='UP')
      CALL mma_allocate(RAW,NRAW,Label='RAW')
      CALL MKRAW(NVERT,DOWN,UP,RAW)
!
!     COMPUTE MIDLEVEL AND LIMITS ON MIDVERTICES
!
      NLTV=NLEV+2
      CALL mma_allocate(LTV,NLTV,Label='LTV')
      CALL MKMID(NVERT,NLEV,DRT,DAW,RAW,LTV,MIDLEV, NMIDV, MIDV1, MIDV2, MXUP, MXDWN)
      CALL mma_deallocate(LTV)
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
      CALL MKCOT(NSYM,NLEV,NVERT,MIDLEV,NMIDV,MIDV1,MIDV2,NWALK,NIPWLK,NSM,DOWN,NOW1,IOW1,NCSF,IOCSF,NOCSF,SCR)
!
!     CONSTRUCT THE CASE LIST
!
      NICASE=NWALK*NIPWLK
      CALL mma_allocate(ICASE,NICASE,Label='ICASE')
      Call MKCLIST(NSYM,NLEV,NVERT,MIDLEV,MIDV1,MIDV2,NMIDV,NICASE,NIPWLK,NSM,DOWN,NOW1,IOW1,ICASE,SCR)
      CALL mma_deallocate(SCR)
!
!     SET UP ENUMERATION TABLES
!
      NUSGN=MXUP*NMIDV
      NLSGN=MXDWN*NMIDV
      CALL mma_allocate(USGN,NUSGN,Label='USGN')
      CALL mma_allocate(LSGN,NLSGN,Label='LSGN')
      Call MKSGNUM(STSYM,NSYM,NLEV,NVERT,MIDLEV,NMIDV,MXUP,MXDWN,NICASE,NIPWLK,DOWN,UP,DAW,RAW,NOW1,IOW1,USGN,LSGN,ICASE)
!
!     EXIT
!
      DRTP => Null()
      DOWNP=> Null()

      END SUBROUTINE MKGUGA

      SUBROUTINE MKGUGA_FREE()
!
!     PURPOSE: FREE THE GUGA TABLES
!
      use stdalloc, only: mma_deallocate
      use gugx, only:  DRT,  DOWN,  UP,  RAW,  DAW,  NOCSF,   &
     &                 IOCSF,  ICASE, USGN, LSGN, NOW1, IOW1, &
     &                 LTV
      IMPLICIT None
!
      If (Allocated(DRT)) Call mma_deallocate(DRT)
      If (Allocated(DOWN)) Call mma_deallocate(DOWN)

      If (Allocated(DAW)) Call mma_deallocate(DAW)
      If (Allocated(UP)) Call mma_deallocate(UP)
      If (Allocated(RAW)) Call mma_deallocate(RAW)
      If (Allocated(LTV)) Call mma_deallocate(LTV)

      If (Allocated(NOW1)) Call mma_deallocate(NOW1)
      If (Allocated(IOW1)) Call mma_deallocate(IOW1)
      If (Allocated(NOCSF)) Call mma_deallocate(NOCSF)
      If (Allocated(IOCSF)) Call mma_deallocate(IOCSF)

      If (Allocated(ICASE)) Call mma_deallocate(ICASE)

      If (Allocated(USGN)) Call mma_deallocate(USGN)
      If (Allocated(LSGN)) Call mma_deallocate(LSGN)

      END SUBROUTINE MKGUGA_FREE
#endif
      END SUBROUTINE GUGACTL
