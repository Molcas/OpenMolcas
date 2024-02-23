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
      use stdalloc, only: mma_allocate, mma_deallocate
      use gugx, only:IA0, IB0, IC0, ISM, LM1RAS, LM3RAS, LV1RAS, LV3RAS,&
     &               MXEO, NVERT,                  NLEV, IFCAS,         &
     &               NCSF, NWALK, MidLev,NMIDV,MVSta,MVEnd,             &
     &                VTAB, DOWN,  ICOUP, IOCP, UP,    MVR, MVL,        &
     &               NVTAB,       NICOUP,NIOCP,       NMVR,NMVL,        &
     &                LTV, MAW, RAW, DAW, NOW1, DRT,  NOCP, NOCSF,      &
     &                    NMAW,                      NNOCP
      IMPLICIT None
#include "rasdim.fh"
#include "caspt2.fh"
#include "segtab.fh"
      Integer, Allocatable:: IVR(:), ISGM(:), NRL(:), ILNDW(:), SCR(:)
      Integer                                NNRL,   NILNDW,   NSCR
      Real*8, Allocatable:: VSGM(:), VTAB_TMP(:), VAL(:)
      Integer                       NVTAB_TMP

      Interface
      SUBROUTINE MKGUGA(NSM,NLEV,NSYM,STSYM,NCSF,Skip_MKSGNUM)
      IMPLICIT None
      Integer NLEV, NSYM, STSYM
      Integer NSM(NLEV)
      Integer NCSF(NSYM)
      Logical, Optional:: Skip_MKSGNUM
      End SUBROUTINE MKGUGA
      SUBROUTINE MKMAW_CP2(IDOWN,IDAW,IUP,IRAW,IMAW,NVERT, MVSta, MVEnd)
      IMPLICIT None
      Integer NVERT, MVSta, MVEnd
      Integer IDOWN(NVERT,0:3),IDAW(NVERT,0:4)
      Integer IUP(NVERT,0:3),IRAW(NVERT,0:4)
      Integer IMAW(NVERT,0:3)
      END SUBROUTINE MKMAW_CP2
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
      IF ((2*IA0+IB0).NE.NACTEL) GOTO 9001
      IF((IA0.LT.0).OR.(IB0.LT.0).OR.(IC0.LT.0)) GOTO 9001

      CALL MKGUGA(ISM,NLEV,NSYM,STSYM,NCSF,Skip_MKSGNUM=.TRUE.)

! DECIDE MIDLEV AND CALCULATE MODIFIED ARC WEIGHT TABLE.

      NMAW=4*NVERT
      CALL mma_allocate(MAW,NMAW,Label='MAW')
      CALL MKMAW_CP2(DOWN,DAW,UP,RAW,MAW,NVERT, MVSta, MVEnd)
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
      CALL MKSEG_CP2(DRT,DOWN,LTV,IVR,MVL,MVR,ISGM,VSGM)

! NIPWLK: NR OF INTEGERS USED TO PACK EACH UP- OR DOWNWALK.
      MXEO=(NLEV*(NLEV+5))/2
      NNOCP=MXEO*NMIDV*NSYM
      NIOCP=NNOCP
      NNRL=(1+MXEO)*NVERT*NSYM
      CALL mma_allocate(NOCP,NNOCP,Label='NOCP')
      CALL mma_allocate(IOCP,NIOCP,Label='IOCP')
      CALL mma_allocate(NRL,NNRL,Label='NRL')
      CALL NRCOUP_CP2(DRT,ISGM,NOW1,NOCP,IOCP,NOCSF,NRL,MVL,MVR)
      CALL mma_deallocate(NRL)

      NILNDW=NWALK
      NVTAB_TMP=20000
      NSCR=7*(NLEV+1)
      CALL mma_allocate(ICOUP,3*NICOUP,Label='ICOUP')
      CALL mma_allocate(VTAB_TMP,NVTAB_TMP,Label='VTAB_TMP')
      CALL mma_allocate(ILNDW,NILNDW,Label='ILNDW')
      CALL mma_allocate(SCR,NSCR,Label='SCR')
      CALL mma_allocate(VAL,NLEV+1,Label='VAL')
      CALL MKCOUP_CP2(nLev,ISm,nVert,MidLev,nMidV,MVSta,MVEnd,          &
     &                IVR,MAW,ISGM,VSGM,NOW1,     NOCP,IOCP,ILNDW,      &
     &                ICOUP,NVTAB_TMP,VTAB_TMP,NVTAB,SCR,VAL)

      CALL mma_allocate(VTAB,NVTAB,Label='VTAB')
      VTAB(1:NVTAB)=VTAB_TMP(1:NVTAB)

      Call mma_deallocate(VTAB_TMP)
      Call mma_deallocate(ILNDW)
      Call mma_deallocate(SCR)
      Call mma_deallocate(VAL)
      Call mma_deallocate(ISGM)
      Call mma_deallocate(VSGM)
      Call mma_deallocate(MAW)
      Call mma_deallocate(IVR)

      CALL mma_deallocate(DRT)
      Call mma_deallocate(DOWN)
      CALL mma_deallocate(DAW)
      CALL mma_deallocate(UP)
      CALL mma_deallocate(RAW)
      Call mma_deallocate(LTV)

      RETURN
 9001 WRITE(6,*)' ERROR IN SUBROUTINE GINIT.'
      WRITE(6,*)'  NR OF ACTIVE ORBITALS:',NLEV
      WRITE(6,*)' NR OF ACTIVE ELECTRONS:',NACTEL
      WRITE(6,*)'        SPIN DEGENERACY:',ISPIN
      CALL ABEND()

      END SUBROUTINE GINIT_CP2

