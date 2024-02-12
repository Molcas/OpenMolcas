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
      SUBROUTINE MKGUGA(NSM,IPRINT)
C
C     PURPOSE: MAKE THE GUGA TABLES
C     NOTE:    TO RETAIN THE TABLES AVAILABLE FOR LATER PURPOSES
C              THE START ADRESSES OF OF THE ARRAYS ARE STORED IN
C              THE COMMON /GUGA/. THESE ARE:
C              DRT,DOWN,DAW,UP,RAW,NOW1,IOW1,NOCSF,IOCSF
C
      use stdalloc, only: mma_allocate, mma_deallocate
      use gugx, only: NLEV, IA0, IB0, IC0, NVERT0,
     &                IFCAS, NVERT, NDRT,  DRT,
     &                NDOWN,  DOWN,  UP, NUP,  RAW, NRAW, MIDLEV,
     &                NMIDV, MXUP, MXDWN, NWALK, NNOW,  DAW, NDAW,
     &                NIOW, NIPWLK, NICASE,  ICASE,       NNOCSF,
     &                 NOCSF, NIOCSF,  IOCSF, LLSGN, LUSGN, NOW1,
     &                IOW1

      IMPLICIT REAL*8 (A-H,O-Z)
C
#include "rasdim.fh"
#include "general.fh"
#include "output_ras.fh"
#include "WrkSpc.fh"
C
      DIMENSION NSM(*)
      Integer, Pointer:: DRTP(:)=>Null(), DOWNP(:)=>Null()
      Integer, Allocatable, Target:: DRT0(:), DOWN0(:)
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
      CALL GETMEM('LTMP','ALLO','INTEGER',LTMP,NTMP)
      CALL mkDRT0 (IA0,IB0,IC0,NVERT0,DRTP,DOWNP,
     *           NTMP,IWORK(LTMP))
      CALL GETMEM('LTMP','FREE','INTEGER',LTMP,NTMP)
C
      IF(IPRINT.ge.DEBUG) THEN
        Write(LF,*)
        Write(LF,*)' PALDUS DRT TABLE (UNRESTRICTED):'
        CALL PRDRT(NVERT0,DRTP,DOWNP)
      ENDIF
C
C     IF THIS IS A RAS CALCULATION PUT UP RESTRICTIONS BY DELETING
C     VERTICES WHICH VIOLATE THE FORMER.
C
      IF(IFCAS.NE.0) THEN
        CALL GETMEM('LV11','ALLO','INTEG',LV,NVERT0)
        CALL RESTR(DRT0,DOWN0,IWORK(LV))
C
C     REASSEMBLE THE DRT TABLE (REMOVE DISCONNECTED VERTICES)
C
        NDRT=5*NVERT
        NDOWN=4*NVERT
        CALL mma_allocate(DRT,NDRT,Label='DRT')
        CALL mma_allocate(DOWN,NDOWN,Label='DOWN')
        CALL mkDRT(DRT0,DOWN0,IWORK(LV),
     *             DRT,DOWN)
        CALL GETMEM('LV11','FREE','INTEG',LV,NVERT0)
        CALL mma_deallocate(DRT0)
        CALL mma_deallocate(DOWN0)
C
        IF(IPRINT.ge.DEBUG) THEN
          Write(LF,*)
          Write(LF,*)' PALDUS DRT TABLE (RESTRICTED):'
          CALL PRDRT(NVERT,DRT,DOWN)
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
      CALL mma_allocate(DAW,NDAW,Label='DAW')
      CALL MKDAW(DOWN,DAW,IPRINT)
C
C     COMPUTE UPCHAIN TABLE AND REVERSE ARC WEIGHTS
C
      NUP=4*NVERT
      NRAW=5*NVERT
      CALL mma_allocate(UP,NUP,Label='UP')
      CALL mma_allocate(RAW,NRAW,Label='RAW')
      CALL MKRAW(DOWN,UP,RAW,IPRINT)
C
C     COMPUTE MIDLEVEL AND LIMITS ON MIDVERTICES
C
      NLTV=NLEV+2
      CALL GETMEM('LTV1','ALLO','INTEG',LLTV,NLTV)
      CALL MKMID(DRT,DAW,RAW,IWORK(LLTV),IPRINT)
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
      CALL mma_allocate(NOCSF,NNOCSF,Label='NOCSF')
      CALL mma_allocate(IOCSF,NIOCSF,Label='IOCSF')
      CALL GETMEM('SCR1','ALLO','INTEG',LSCR,NSCR)
      CALL MKCOT(NSM,DOWN,NOW1,IOW1,
     *           IOCSF,NOCSF,IWORK(LSCR),IPRINT)
C
C     CONSTRUCT THE CASE LIST
C
      NICASE=NWALK*NIPWLK
      CALL mma_allocate(ICASE,NICASE,Label='ICASE')
      CALL MKCLIST(NSM,DOWN,NOW1,IOW1,
     &            ICASE,IWORK(LSCR))
      CALL GETMEM('SCR1','FREE','INTEG',LSCR,NSCR)
C
C     SET UP ENUMERATION TABLES
C
      NUSGN=MXUP*NMIDV
      NLSGN=MXDWN*NMIDV
      CALL GETMEM('IUSG','ALLO','INTEG',LUSGN,NUSGN)
      CALL GETMEM('ILSG','ALLO','INTEG',LLSGN,NLSGN)
      CALL MKSGNUM(DOWN,UP,DAW,RAW,
     *             NOW1,IOW1,IWORK(LUSGN),IWORK(LLSGN),
     *             ICASE,IPRINT)
C
C     EXIT
C
      DRTP => Null()
      DOWNP=> Null()
      END

      SUBROUTINE MKGUGA_FREE()
C
C     PURPOSE: FREE THE GUGA TABLES
C
      use stdalloc, only: mma_deallocate
      use gugx, only:  DRT,  DOWN,  UP,  RAW,  DAW,  NOCSF,
     &                 IOCSF,  ICASE, LUSGN, LLSGN, NOW1, IOW1,
     &                MXUP, MXDWN, NMIDV
      IMPLICIT REAL*8 (A-H,O-Z)
C

      Call mma_deallocate(DRT)
      Call mma_deallocate(DOWN)

      Call mma_deallocate(DAW)
      Call mma_deallocate(UP)
      Call mma_deallocate(RAW)

      Call mma_deallocate(NOW1)
      Call mma_deallocate(IOW1)
      Call mma_deallocate(NOCSF)
      Call mma_deallocate(IOCSF)

      Call mma_deallocate(ICASE)

      NUSGN=MXUP*NMIDV
      NLSGN=MXDWN*NMIDV
      CALL GETMEM('IUSG','FREE','INTEG',LUSGN,NUSGN)
      CALL GETMEM('ILSG','FREE','INTEG',LLSGN,NLSGN)

      END
