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
      IMPLICIT REAL*8 (A-H,O-Z)
#include "rasdim.fh"
#include "caspt2.fh"
#include "WrkSpc.fh"
#include "pt2_guga.fh"
#include "segtab.fh"
      Integer, Allocatable, Target:: DRT0(:), DRT(:)
      Integer, Allocatable, Target:: DOWN0(:), DOWN(:)
      Integer, Pointer:: DRTP(:)=>Null(), DOWNP(:)=>Null()
      Integer, Allocatable:: TMP(:), V11(:), DAW(:), LTV(:)


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
      CALL GETMEM('UP','ALLO','INTEG',LUP,NUP)
C REVERSE ARC WEIGHT TABLE:
      NRAW=5*NVERT
      CALL GETMEM('RAW','ALLO','INTEG',LRAW,NRAW)
C DECIDE MIDLEV AND CALCULATE MODIFIED ARC WEIGHT TABLE.
      NMAW=4*NVERT
      CALL GETMEM('MAW','ALLO','INTEG',LMAW,NMAW)
      CALL MKMAW_CP2(DOWN,DAW,IWORK(LUP),IWORK(LRAW),
     &           IWORK(LMAW),LTV)
C THE DAW, UP AND RAW TABLES WILL NOT BE NEEDED ANY MORE:
      CALL mma_deallocate(DAW)
      CALL GETMEM('UP','FREE','INTEG',LUP,NUP)
      CALL GETMEM('RAW','FREE','INTEG',LRAW,NRAW)
C CALCULATE SEGMENT VALUES. ALSO, MVL AND MVR TABLES.
      NIVR=2*NVERT
      CALL GETMEM('IVR','ALLO','INTEG',LIVR,NIVR)
      NMVL=2*NMIDV
      NMVR=2*NMIDV
      CALL GETMEM('MVL','ALLO','INTEG',LMVL,NMVL)
      CALL GETMEM('MVR','ALLO','INTEG',LMVR,NMVR)
      NSGMNT=26*NVERT
      CALL GETMEM('ISGM','ALLO','INTEG',LISGM,NSGMNT)
      CALL GETMEM('VSGM','ALLO','REAL',LVSGM,NSGMNT)
      CALL MKSEG_CP2(DRT,DOWN,LTV,IWORK(LIVR),
     &           IWORK(LMVL),IWORK(LMVR),IWORK(LISGM),WORK(LVSGM))
      Call mma_deallocate(DOWN)
      Call mma_deallocate(LTV)
C FORM VARIOUS OFFSET TABLES:
      NNOW=2*NMIDV*NSYM
      NIOW=NNOW
      CALL GETMEM('NOW','ALLO','INTEG',LNOW,NNOW)
      CALL GETMEM('IOW','ALLO','INTEG',LIOW,NIOW)
      MXEO=(NLEV*(NLEV+5))/2
      NNOCP=MXEO*NMIDV*NSYM
      NIOCP=NNOCP
      NNRL=(1+MXEO)*NVERT*NSYM
      CALL GETMEM('NOCP','ALLO','INTEG',LNOCP,NNOCP)
      CALL GETMEM('IOCP','ALLO','INTEG',LIOCP,NIOCP)
      CALL GETMEM('NRL','ALLO','INTEG',LNRL,NNRL)
      NNOCSF=NMIDV*(NSYM**2)
      NIOCSF=NNOCSF
C NIPWLK: NR OF INTEGERS USED TO PACK EACH UP- OR DOWNWALK.
      NIPWLK=1+(MIDLEV-1)/15
      NIPWLK=MAX(NIPWLK,1+(NLEV-MIDLEV-1)/15)
      CALL GETMEM('NOCSF','ALLO','INTEG',LNOCSF,NNOCSF)
      CALL GETMEM('IOCSF','ALLO','INTEG',LIOCSF,NIOCSF)
      CALL NRCOUP_CP2(DRT,IWORK(LISGM),IWORK(LNOW),
     &            IWORK(LIOW),IWORK(LNOCP),IWORK(LIOCP),IWORK(LNOCSF),
     &            IWORK(LIOCSF),IWORK(LNRL),IWORK(LMVL),IWORK(LMVR))
      CALL mma_deallocate(DRT)
      CALL GETMEM('NRL','FREE','INTEG',LNRL,NNRL)
      NILNDW=NWALK
      NICASE=NWALK*NIPWLK
      CALL GETMEM('ICASE','ALLO','INTEG',LICASE,NICASE)
      NNICOUP=3*NICOUP
      CALL GETMEM('ICOUP','ALLO','INTEG',LICOUP,NNICOUP)
      NVTAB_TMP=20000
      CALL GETMEM('VTAB_TMP','ALLO','REAL',LVTAB_TMP,NVTAB_TMP)
      NSCR=7*(NLEV+1)
      CALL GETMEM('ILNDW','ALLO','INTEG',LILNDW,NILNDW)
      CALL GETMEM('SCR','ALLO','INTEG',LSCR,NSCR)
      CALL GETMEM('VAL','ALLO','REAL',LVAL,NLEV+1)
      CALL MKCOUP_CP2(IWORK(LIVR),IWORK(LMAW),IWORK(LISGM),
     &            WORK(LVSGM),IWORK(LNOW),IWORK(LIOW),IWORK(LNOCP),
     &     IWORK(LIOCP),IWORK(LILNDW),IWORK(LICASE),IWORK(LICOUP),
     &     NVTAB_TMP,WORK(LVTAB_TMP),NVTAB_FINAL,IWORK(LSCR),
     &     WORK(LVAL))
* Set NVTAB in common block /IGUGA/ in file pt2_guga.fh:
      NVTAB=NVTAB_FINAL
      CALL GETMEM('VTAB','ALLO','REAL',LVTAB,NVTAB)
      CALL DCOPY_(NVTAB,WORK(LVTAB_TMP),1,WORK(LVTAB),1)
      CALL GETMEM('VTAB_TMP','FREE','REAL',LVTAB_TMP,NVTAB_TMP)
      CALL GETMEM('ILNDW','FREE','INTEG',LILNDW,NILNDW)
      CALL GETMEM('SCR','FREE','INTEG',LSCR,NSCR)
      CALL GETMEM('VAL','FREE','REAL',LVAL,NLEV+1)
      CALL GETMEM('ISGM','FREE','INTEG',LISGM,NSGMNT)
      CALL GETMEM('VSGM','FREE','REAL',LVSGM,NSGMNT)
      CALL GETMEM('MAW','FREE','INTEG',LMAW,NMAW)
      CALL GETMEM('IVR','FREE','INTEG',LIVR,NIVR)

      DOWNP=>Null()
      DRTP=>Null()

      RETURN
 9001 WRITE(6,*)' ERROR IN SUBROUTINE GINIT.'
      WRITE(6,*)'  NR OF ACTIVE ORBITALS:',NLEV
      WRITE(6,*)' NR OF ACTIVE ELECTRONS:',NACTEL
      WRITE(6,*)'        SPIN DEGENERACY:',ISPIN
      CALL ABEND()
      END

