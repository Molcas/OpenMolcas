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
      SUBROUTINE MKGUGA_m(NSM,IPRINT)
C
C     PURPOSE: MAKE THE GUGA TABLES
C     NOTE:    TO RETAIN THE TABLES AVAILABLE FOR LATER PURPOSES
C              THE START ADRESSES OF OF THE ARRAYS ARE STORED IN
C              THE COMMON /GUGA/. THESE ARE:
C              DRT,LDOWN,LDAW,LUP,LRAW,NOW1,IOW1,LNOCSF,LIOCSF
C
      use mcpdft_output, only:  debug, lf
      use stdalloc, only: mma_allocate, mma_deallocate
      use gugx, only: NLEV, IA0, IB0, IC0, NVERT0,
     &                IFCAS, NVERT, NDRT,  DRT,
     &                NDOWN, LDOWN, LUP, NUP, LRAW, NRAW, MIDLEV,
     &                NMIDV, MXUP, MXDWN, NWALK, NNOW, LDAW, NDAW,
     &                NIOW, NIPWLK, NICASE,  ICASE,       NNOCSF,
     &                LNOCSF, NIOCSF, LIOCSF, LLSGN, LUSGN, NOW1,
     &                IOW1

      IMPLICIT REAL*8 (A-H,O-Z)
C
#include "rasdim.fh"
#include "general.fh"
#include "WrkSpc.fh"
C
      DIMENSION NSM(*)
      Integer, Pointer:: DRTP(:)=>Null()
      Integer, Allocatable, Target:: DRT0(:)
C
C     SET UP A FULL PALDUS DRT TABLE:
C     (INITIALLY NO RESTRICTIONS ARE PUT UP)
C
      IAC=MIN(IA0,IC0)
      NVERT0=((IA0+1)*(IC0+1)*(2*IB0+IAC+2))/2-(IAC*(IAC+1)*(IAC+2))/6
      NDRT0=5*NVERT0
      NDOWN0=4*NVERT0
      NTMP=((NLEV+1)*(NLEV+2))/2

      IF(IFCAS.NE.0) THEN
         CALL mma_allocate(DRT0,NDRT0,Label='DRT0')
         CALL GETMEM('DOWN','ALLO','INTEGER',LDOWN0,NDOWN0)
         DRTP => DRT0
         LDOWNP=LDOWN0
      ELSE
         NDRT=NDRT0
         NDOWN=NDOWN0
         CALL mma_allocate(DRT,NDRT,Label='DRT')
         CALL GETMEM('DOWN','ALLO','INTEGER',LDOWN,NDOWN)
         DRTP => DRT
         LDOWNP=LDOWN
         NVERT=NVERT0
      ENDIF

      CALL GETMEM('LTMP','ALLO','INTEGER',LTMP,NTMP)
      CALL mkDRT0(IA0,IB0,IC0,NVERT0,DRTP,IWORK(LDOWNP),
     *           NTMP,IWORK(LTMP))
      CALL GETMEM('LTMP','FREE','INTEGER',LTMP,NTMP)
C
      IF(IPRINT >= DEBUG) THEN
        Write(LF,*)
        Write(LF,*)' PALDUS DRT TABLE (UNRESTRICTED):'
        CALL print_drt(NVERT0,DRTP,IWORK(LDOWNP))
      ENDIF
C
C     IF THIS IS A RAS CALCULATION PUT UP RESTRICTIONS BY DELETING
C     VERTICES WHICH VIOLATE THE FORMER.
C
      IF(IFCAS.NE.0) THEN
        CALL GETMEM('LV11','ALLO','INTEG',LV,NVERT0)
        CALL RESTR_m(DRT0,IWORK(LDOWN0),IWORK(LV))
C
C     REASSEMBLE THE DRT TABLE (REMOVE DISCONNECTED VERTICES)
C
        NDRT=5*NVERT
        NDOWN=4*NVERT
        CALL mma_allocate(DRT,NDRT,Label='DRT')
        CALL GETMEM('DWN1','ALLO','INTEG',LDOWN,NDOWN)
        CALL mkDRT(DRT0,IWORK(LDOWN0),IWORK(LV),
     *             DRT,IWORK(LDOWN))
        CALL GETMEM('LV11','FREE','INTEG',LV,NVERT0)
        CALL mma_deallocate(DRT0)
        CALL GETMEM('DOWN','FREE','INTEG',LDOWN0,NDOWN0)
C
        IF(IPRINT >= DEBUG) THEN
          Write(LF,*)
          Write(LF,*)' PALDUS DRT TABLE (RESTRICTED):'
          CALL print_drt(NVERT,DRT,IWORK(LDOWN))
        ENDIF
C
C     IF THIS IS A CAS CALCULATION PROCEED WITH THE UNRESTRICTED
C     DRT TABLE
C
      ENDIF
C
C     CALCULATE ARC WEIGHT AND LTV TABLES.
C
      NDAW=5*NVERT
      CALL GETMEM('DAW1','ALLO','INTEG',LDAW,NDAW)
      CALL MKDAW_m(IWORK(LDOWN),IWORK(LDAW),IPRINT)
C
C     COMPUTE UPCHAIN TABLE AND REVERSE ARC WEIGHTS
C
      NUP=4*NVERT
      NRAW=5*NVERT
      CALL GETMEM('LUP1','ALLO','INTEG',LUP,NUP)
      CALL GETMEM('RAW1','ALLO','INTEG',LRAW,NRAW)
      CALL MKRAW_m(IWORK(LDOWN),IWORK(LUP),IWORK(LRAW),IPRINT)
C
C     COMPUTE MIDLEVEL AND LIMITS ON MIDVERTICES
C
      NLTV=NLEV+2
      CALL GETMEM('LTV1','ALLO','INTEG',LLTV,NLTV)
      CALL MKMID_m(DRT,IWORK(LDAW),IWORK(LRAW),IWORK(LLTV),
     & IPRINT)
      CALL GETMEM('LTV1','FREE','INTEG',LLTV,NLTV)
C
C     FORM VARIOUS OFFSET TABLES:
C     NOTE: NIPWLK AND DOWNWLK ARE THE NUMER OF INTEGER WORDS USED
C           TO STORE THE UPPER AND LOWER WALKS IN PACKED FORM.
C
      NIPWLK=1+(MIDLEV-1)/15
      NIPWLK=MAX(NIPWLK,1+(NLEV-MIDLEV-1)/15)
      NNOW=2*NMIDV*NSYM
      NIOW=NNOW
      NNOCSF=NMIDV*(NSYM**2)
      NIOCSF=NNOCSF
      NSCR=MAX(6,3*(NLEV+1))
      CALL mma_allocate(NOW1,NNOW,Label='NOW1')
      CALL mma_allocate(IOW1,NIOW,Label='IOW1')
      CALL GETMEM('NCSF','ALLO','INTEG',LNOCSF,NNOCSF)
      CALL GETMEM('ICSF','ALLO','INTEG',LIOCSF,NIOCSF)
      CALL GETMEM('SCR1','ALLO','INTEG',LSCR,NSCR)
      CALL MKCOT_m(NSM,IWORK(LDOWN),NOW1,IOW1,
     *           IWORK(LIOCSF),IWORK(LNOCSF),IWORK(LSCR),IPRINT)
C
C     CONSTRUCT THE CASE LIST
C
      NICASE=NWALK*NIPWLK
      CALL mma_allocate(ICASE,NICASE,Label='ICASE')
      CALL MKCLIST(NSM,IWORK(LDOWN),NOW1,IOW1,
     &             ICASE,IWORK(LSCR))
      CALL GETMEM('SCR1','FREE','INTEG',LSCR,NSCR)
C
C     SET UP ENUMERATION TABLES
C
      NUSGN=MXUP*NMIDV
      NLSGN=MXDWN*NMIDV
      CALL GETMEM('IUSG','ALLO','INTEG',LUSGN,NUSGN)
      CALL GETMEM('ILSG','ALLO','INTEG',LLSGN,NLSGN)
      CALL MKSGNUM_m(IWORK(LDOWN),IWORK(LUP),IWORK(LDAW),IWORK(LRAW),
     *             NOW1,IOW1,IWORK(LUSGN),IWORK(LLSGN),
     *             ICASE,IPRINT)
C
C     EXIT
C
      DRTP => Null()
      END

      SUBROUTINE MKGUGA_FREE_m()
      use stdalloc, only: mma_deallocate
      use gugx, only:  DRT, LDOWN, LUP, LRAW, LDAW, LNOCSF,
     &                LIOCSF,  ICASE, LUSGN, LLSGN, NOW1, IOW1,
     &                MXUP, MXDWN, NMIDV,       NDOWN, NDAW, NUP,
     &                NRAW, NNOCSF, NIOCSF
C
C     PURPOSE: FREE THE GUGA TABLES
C
      IMPLICIT REAL*8 (A-H,O-Z)
C

      CALL mma_deallocate(DRT)
      CALL GETMEM('DWN0/1','FREE','INTE',LDOWN,NDOWN)

      CALL GETMEM('DAW1','FREE','INTE',LDAW,NDAW)
      CALL GETMEM('LUP1','FREE','INTE',LUP,NUP)
      CALL GETMEM('RAW1','FREE','INTE',LRAW,NRAW)

      Call mma_deallocate(NOW1)
      Call mma_deallocate(IOW1)
      CALL GETMEM('NCSF','FREE','INTEG',LNOCSF,NNOCSF)
      CALL GETMEM('ICSF','FREE','INTEG',LIOCSF,NIOCSF)

      Call mma_deallocate(ICASE)

      NUSGN=MXUP*NMIDV
      NLSGN=MXDWN*NMIDV
      CALL GETMEM('IUSG','FREE','INTEG',LUSGN,NUSGN)
      CALL GETMEM('ILSG','FREE','INTEG',LLSGN,NLSGN)

      END
