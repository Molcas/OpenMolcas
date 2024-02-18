************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 1994,2006, Per Ake Malmqvist                           *
************************************************************************
*--------------------------------------------*
* 1994  PER-AAKE MALMQUIST                   *
* DEPARTMENT OF THEORETICAL CHEMISTRY        *
* UNIVERSITY OF LUND                         *
* SWEDEN                                     *
* 2006  PER-AAKE MALMQUIST                   *
*--------------------------------------------*
      SUBROUTINE GINIT_CP2
      use stdalloc, only: mma_allocate, mma_deallocate
      use pt2_guga_data
      IMPLICIT REAL*8 (A-H,O-Z)
#include "rasdim.fh"
#include "caspt2.fh"
#include "segtab.fh"
      Integer, Allocatable, Target:: DRT0(:), DRT(:)
      Integer, Allocatable, Target:: DOWN0(:), DOWN(:)
      Integer, Pointer:: DRTP(:)=>Null(), DOWNP(:)=>Null()
      Integer, Allocatable:: TMP(:), V11(:), DAW(:), LTV(:), RAW(:),
     &                       UP(:), MAW(:), IVR(:), ISGM(:), NRL(:),
     &                       ILNDW(:), SCR(:)
      Real*8, Allocatable:: VSGM(:), VTAB_TMP(:), VAL(:)
      Integer, External:: ip_of_iWork, ip_of_Work

      LV1RAS=NRAS1T
      LV3RAS=LV1RAS+NRAS2T
      LM1RAS=2*LV1RAS-NHOLE1
      LM3RAS=NACTEL-NELE3

C SET UP A FULL PALDUS DRT TABLE:
      IB0=ISPIN-1
      IA0=(NACTEL-IB0)/2
      IC0=NLEV-IA0-IB0
      IF ((2*IA0+IB0).NE.NACTEL) GOTO 9001
      IF((IA0.LT.0).OR.(IB0.LT.0).OR.(IC0.LT.0)) GOTO 9001
      IAC=MIN(IA0,IC0)
      NVERT0=((IA0+1)*(IC0+1)*(2*IB0+IAC+2))/2-(IAC*(IAC+1)*(IAC+2))/6
      NDRT0=5*NVERT0
      NDOWN0=4*NVERT0
      NTMP=((NLEV+1)*(NLEV+2))/2

      IF((NRAS1T+NRAS3T).NE.0) THEN
        CALL mma_allocate(DRT0,NDRT0,Label='DRT0')
        CALL mma_allocate(DOWN0,NDOWN0,Label='DOWN0')
        DRTP=>DRT0
        DOWNP=>DOWN0
      ELSE
        NDRT=NDRT0
        NDOWN=NDOWN0
        CALL mma_allocate(DRT,NDRT,Label='DRT')
        CALL mma_allocate(DOWN,NDOWN,Label='DOWN')
        DRTP=>DRT
        DOWNP=>DOWN
        NVERT=NVERT0
      ENDIF

      CALL mma_allocate(TMP,NTMP,Label='TMP')
      CALL mkDRT0(IA0,IB0,IC0,NVERT0,DRTP,DOWNP,NTMP,TMP)
      CALL mma_deallocate(TMP)
#ifdef _DEBUGPRINT_
      WRITE(6,*)
      WRITE(6,*)' PALDUS DRT TABLE (UNRESTRICTED):'
      CALL PRDRT(NVERT0,DRTP,DOWNP)
#endif
C RESTRICTIONS?
      IF((NRAS1T+NRAS3T).NE.0) THEN
C  CONSTRUCT A RESTRICTED GRAPH.
        CALL mma_allocate(V11,NVERT0,Label='V11')
        CALL RESTR(NVERT0,DRTP,DOWN0,V11,
     &                 LV1RAS,LV3RAS,LM1RAS,LM3RAS,NVERT)

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
        WRITE(6,*)
        WRITE(6,*)' PALDUS DRT TABLE (RESTRICTED):'
        CALL PRDRT(NVERT,DRT,DOWN)
#endif
      END IF

C CALCULATE DIRECT ARC WEIGHT AND LTV TABLES.
      NDAW=5*NVERT
      CALL mma_allocate(DAW,NDAW,Label='DAW')
      NLTV=NLEV+2
      CALL mma_allocate(LTV,NLTV,Label='LTV')
      CALL MKDAW_CP2(NLEV,NVERT,DRT,DOWN,DAW,LTV)
C UPCHAIN INDEX TABLE:
      NUP=4*NVERT
      CALL mma_allocate(UP,NUP,Label='UP')
C REVERSE ARC WEIGHT TABLE:
      NRAW=5*NVERT
      CALL mma_allocate(RAW,NRAW,Label='RAW')
C DECIDE MIDLEV AND CALCULATE MODIFIED ARC WEIGHT TABLE.
      NMAW=4*NVERT
      CALL mma_allocate(MAW,NMAW,Label='MAW')
      CALL MKMAW_CP2(DOWN,DAW,UP,RAW,MAW,LTV)
C THE DAW, UP AND RAW TABLES WILL NOT BE NEEDED ANY MORE:
      CALL mma_deallocate(DAW)
      CALL mma_deallocate(UP)
      CALL mma_deallocate(RAW)
C CALCULATE SEGMENT VALUES. ALSO, MVL AND MVR TABLES.
      NIVR=2*NVERT
      CALL mma_allocate(IVR,NIVR,Label='IVR')
      NMVL=2*NMIDV
      NMVR=2*NMIDV
      CALL mma_allocate(MVL,NMVL,Label='MVL')
      CALL mma_allocate(MVR,NMVR,Label='MVR')
      NSGMNT=26*NVERT
      CALL mma_allocate(ISGM,NSGMNT,Label='ISGM')
      CALL mma_allocate(VSGM,NSGMNT,Label='VSGM')
      CALL MKSEG_CP2(DRT,DOWN,LTV,IVR,MVL,MVR,ISGM,VSGM)
      Call mma_deallocate(DOWN)
      Call mma_deallocate(LTV)
C FORM VARIOUS OFFSET TABLES:
      NNOW=2*NMIDV*NSYM
      NIOW=NNOW
      CALL mma_allocate(NOW1,NNOW,Label='NOW1')
      LNOW = ip_of_iWork(NOW1(1))
      CALL mma_allocate(IOW1,NIOW,Label='IOW1')
      LIOW = ip_of_iWork(IOW1(1))
      MXEO=(NLEV*(NLEV+5))/2
      NNOCP=MXEO*NMIDV*NSYM
      NIOCP=NNOCP
      NNRL=(1+MXEO)*NVERT*NSYM
      CALL mma_allocate(NOCP,NNOCP,Label='NOCP')
      LNOCP = ip_of_iWork(NOCP(1))
      CALL mma_allocate(IOCP,NIOCP,Label='IOCP')
      LIOCP = ip_of_iWork(IOCP(1))
      CALL mma_allocate(NRL,NNRL,Label='NRL')
      NNOCSF=NMIDV*(NSYM**2)
      NIOCSF=NNOCSF
C NIPWLK: NR OF INTEGERS USED TO PACK EACH UP- OR DOWNWALK.
      NIPWLK=1+(MIDLEV-1)/15
      NIPWLK=MAX(NIPWLK,1+(NLEV-MIDLEV-1)/15)
      CALL mma_allocate(NOCSF,NNOCSF,Label='NOCSF')
      LNOCSF = ip_of_iWork(NOCSF(1))
      CALL mma_allocate(IOCSF,NIOCSF,Label='IOCSF')
      LIOCSF = ip_of_iWork(IOCSF(1))
      CALL NRCOUP_CP2(DRT,ISGM,NOW1,IOW1,NOCP,IOCP,NOCSF,IOCSF,
     &                NRL,MVL,MVR)
      CALL mma_deallocate(DRT)
      CALL mma_deallocate(NRL)
      NILNDW=NWALK
      NICASE=NWALK*NIPWLK
      CALL mma_allocate(ICASE,NICASE,Label='ICASE')
      LICASE = ip_of_iWork(ICASE(1))
      NNICOUP=3*NICOUP
      CALL mma_allocate(ICOUP,NNICOUP,Label='ICOUP')
      LICOUP = ip_of_iWork(ICOUP(1))
      NVTAB_TMP=20000
      CALL mma_allocate(VTAB_TMP,NVTAB_TMP,Label='VTAB_TMP')
      NSCR=7*(NLEV+1)
      CALL mma_allocate(ILNDW,NILNDW,Label='ILNDW')
      CALL mma_allocate(SCR,NSCR,Label='SCR')
      CALL mma_allocate(VAL,NLEV+1,Label='VAL')
      CALL MKCOUP_CP2(IVR,MAW,ISGM,VSGM,NOW1,IOW1,NOCP,IOCP,ILNDW,ICASE,
     &                ICOUP,NVTAB_TMP,VTAB_TMP,NVTAB_FINAL,SCR,VAL)

* Set NVTAB in common block /IGUGA/ in file pt2_guga.fh:
      NVTAB=NVTAB_FINAL
      CALL mma_allocate(VTAB,NVTAB,Label='VTAB')
      LVTAB = ip_of_Work(VTAB(1))
      VTAB(1:NVTAB)=VTAB_TMP(1:NVTAB)

      Call mma_deallocate(VTAB_TMP)
      Call mma_deallocate(ILNDW)
      Call mma_deallocate(SCR)
      Call mma_deallocate(VAL)
      Call mma_deallocate(ISGM)
      Call mma_deallocate(VSGM)
      Call mma_deallocate(MAW)
      Call mma_deallocate(IVR)

      DOWNP=>Null()
      DRTP=>Null()

      RETURN
 9001 WRITE(6,*)' ERROR IN SUBROUTINE GINIT.'
      WRITE(6,*)'  NR OF ACTIVE ORBITALS:',NLEV
      WRITE(6,*)' NR OF ACTIVE ELECTRONS:',NACTEL
      WRITE(6,*)'        SPIN DEGENERACY:',ISPIN
      CALL ABEND()
      END

