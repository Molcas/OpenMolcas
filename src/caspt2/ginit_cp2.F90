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
! Copyright (C) 1994,2006, Per Ake Malmqvist                           *
!***********************************************************************
!--------------------------------------------*
! 1994  PER-AAKE MALMQUIST                   *
! DEPARTMENT OF THEORETICAL CHEMISTRY        *
! UNIVERSITY OF LUND                         *
! SWEDEN                                     *
! 2006  PER-AAKE MALMQUIST                   *
!--------------------------------------------*
      SUBROUTINE GINIT_CP2()
      use Definitions, only: u6
      use stdalloc, only: mma_allocate, mma_deallocate
      use gugx, only:IA0, IB0, IC0, ISM, LM1RAS, LM3RAS, LV1RAS, LV3RAS,&
     &               MXEO,                         NLEV, IFCAS,         &
     &               NCSF, NWALK, NMIDV,                         &
     &                VTAB,        ICOUP, IOCP,        MVR, MVL,        &
     &               NVTAB,       NICOUP,NIOCP,       NMVR,NMVL,        &
     &                          RAW, DAW, NOW1,       NOCP, NOCSF,      &
     &                                    IOW1,      NNOCP,             &
     &               ICASE, NICASE, SGS
      IMPLICIT None
#include "rasdim.fh"
#include "caspt2.fh"
#include "segtab.fh"
      Integer, Allocatable:: IVR(:), ISGM(:), NRL(:), ILNDW(:), SCR(:)
      Integer                                NNRL,   NILNDW,   NSCR
      Real*8, Allocatable:: VSGM(:), VTAB_TMP(:), VAL(:)
      Integer                       NVTAB_TMP, NVERT

      Interface
      SUBROUTINE MKGUGA(NLEV,NSYM,STSYM,NCSF,Skip_MKSGNUM)
      IMPLICIT None
      Integer NLEV, NSYM, STSYM
      Integer NCSF(NSYM)
      Logical, Optional:: Skip_MKSGNUM
      End SUBROUTINE MKGUGA
      End Interface

      LV1RAS=NRAS1T
      LV3RAS=LV1RAS+NRAS2T
      LM1RAS=2*LV1RAS-NHOLE1
      LM3RAS=NACTEL-NELE3
      IF ((NRAS1T+NRAS3T)/=0) Then
         IFCAS=1
      Else
         IFCAS=0
      End If

! SET UP A FULL PALDUS DRT TABLE:
      IB0=ISPIN-1
      IA0=(NACTEL-IB0)/2
      IC0=NLEV-IA0-IB0
      IF ( ((2*IA0+IB0).NE.NACTEL) .or.                    &
         ((IA0.LT.0).OR.(IB0.LT.0).OR.(IC0.LT.0)) ) Then
         WRITE(u6,*)' ERROR IN SUBROUTINE GINIT.'
         WRITE(u6,*)'  NR OF ACTIVE ORBITALS:',NLEV
         WRITE(u6,*)' NR OF ACTIVE ELECTRONS:',NACTEL
         WRITE(u6,*)'        SPIN DEGENERACY:',ISPIN
        CALL ABEND()
      End If

      CALL MKGUGA(NLEV,NSYM,STSYM,NCSF,Skip_MKSGNUM=.TRUE.)
      nVert=SGS%nVert

! DECIDE MIDLEV AND CALCULATE MODIFIED ARC WEIGHT TABLE.

      CALL mma_allocate(SGS%MAW,4*NVERT,Label='MAW')
      CALL MKMAW(SGS%DOWN,DAW,SGS%UP,RAW,SGS%MAW,NVERT, SGS%MVSta, SGS%MVEnd)
! THE DAW, UP AND RAW TABLES WILL NOT BE NEEDED ANY MORE:

! CALCULATE SEGMENT VALUES. ALSO, MVL AND MVR TABLES.
      NIVR=2*NVERT
      NMVL=2*NMIDV
      NMVR=2*NMIDV
      NSGMNT=26*NVERT
      CALL mma_allocate(IVR,NIVR,Label='IVR')
      CALL mma_allocate(MVL,NMVL,Label='MVL')
      CALL mma_allocate(MVR,NMVR,Label='MVR')
      CALL mma_allocate(ISGM,NSGMNT,Label='ISGM')
      CALL mma_allocate(VSGM,NSGMNT,Label='VSGM')
      CALL MKSEG_CP2(SGS%DRT,SGS%DOWN,SGS%LTV,IVR,MVL,MVR,ISGM,VSGM,nVert)

! NIPWLK: NR OF INTEGERS USED TO PACK EACH UP- OR DOWNWALK.
      MXEO=(NLEV*(NLEV+5))/2
      NNOCP=MXEO*NMIDV*NSYM
      NIOCP=NNOCP
      NNRL=(1+MXEO)*NVERT*NSYM
      CALL mma_allocate(NOCP,NNOCP,Label='NOCP')
      CALL mma_allocate(IOCP,NIOCP,Label='IOCP')
      CALL mma_allocate(NRL,NNRL,Label='NRL')
      CALL NRCOUP_CP2(SGS%DRT,ISGM,NOW1,NOCP,IOCP,NOCSF,NRL,MVL,MVR,nVert)
      CALL mma_deallocate(NRL)

      NILNDW=NWALK
      NVTAB_TMP=20000
      NSCR=7*(NLEV+1)
      CALL mma_allocate(ICOUP,3*NICOUP,Label='ICOUP')
      CALL mma_allocate(VTAB_TMP,NVTAB_TMP,Label='VTAB_TMP')
      CALL mma_allocate(ILNDW,NILNDW,Label='ILNDW')
      CALL mma_allocate(SCR,NSCR,Label='SCR')
      CALL mma_allocate(VAL,NLEV+1,Label='VAL')
      CALL MKCOUP(nSym,nLev,ISm,nVert,SGS%MidLev,nMidV,SGS%MVSta,SGS%MVEnd,     &
     &            MxEO,nICoup,nWalk,nICase,nVTAB_TMP,               &
     &            IVR,SGS%MAW,ISGM,VSGM,NOW1,IOW1,NOCP,IOCP,ILNDW,      &
     &            ICase, ICOUP,VTAB_TMP,NVTAB,SCR,VAL)

      CALL mma_allocate(VTAB,NVTAB,Label='VTAB')
      VTAB(1:NVTAB)=VTAB_TMP(1:NVTAB)

      Call mma_deallocate(VTAB_TMP)
      Call mma_deallocate(ILNDW)
      Call mma_deallocate(SCR)
      Call mma_deallocate(VAL)
      Call mma_deallocate(ISGM)
      Call mma_deallocate(VSGM)
      Call mma_deallocate(SGS%MAW)
      Call mma_deallocate(IVR)

      CALL mma_deallocate(SGS%DRT)
      Call mma_deallocate(SGS%DOWN)
      CALL mma_deallocate(DAW)
      CALL mma_deallocate(SGS%UP)
      CALL mma_deallocate(RAW)
      Call mma_deallocate(SGS%LTV)

      RETURN

      END SUBROUTINE GINIT_CP2

