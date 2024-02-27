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
      use gugx, only:IA0, IB0, IC0,      LM1RAS, LM3RAS, LV1RAS, LV3RAS,&
     &               IFCAS, CIS, SGS, EXS
      IMPLICIT None
#include "rasdim.fh"
#include "caspt2.fh"
#include "segtab.fh"
      Integer, Allocatable:: IVR(:), ISGM(:), NRL(:), ILNDW(:), SCR(:)
      Integer                                NNRL,   NILNDW,   NSCR
      Real*8, Allocatable:: VSGM(:), VTAB_TMP(:), VAL(:)
      Integer NICOUP, NVTAB_TMP, NICASE, nVTab

      Interface
      SUBROUTINE MKGUGA(NLEV,NSYM,STSYM,Skip_MKSGNUM)
      IMPLICIT None
      Integer NLEV, NSYM, STSYM
      Logical, Optional:: Skip_MKSGNUM
      End SUBROUTINE MKGUGA
      End Interface

      Associate ( nLev => SGS%nLev, nWalk => CIS%nWalk,                 &
     &            nVert=> SGS%nVert, nMidV=>CIS%nMidV, MXEO => EXS%MxEO)

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

      CALL MKGUGA(NLEV,NSYM,STSYM,Skip_MKSGNUM=.TRUE.)


! DECIDE MIDLEV AND CALCULATE MODIFIED ARC WEIGHT TABLE.

      CALL mma_allocate(SGS%MAW,4*NVERT,Label='MAW')
      CALL MKMAW(SGS%DOWN,SGS%DAW,SGS%UP,SGS%RAW,SGS%MAW,NVERT, SGS%MVSta, SGS%MVEnd)
! THE DAW, UP AND RAW TABLES WILL NOT BE NEEDED ANY MORE:

! CALCULATE SEGMENT VALUES. ALSO, MVL AND MVR TABLES.
      NIVR=2*NVERT
      NSGMNT=26*NVERT
      CALL mma_allocate(IVR,NIVR,Label='IVR')
      CALL mma_allocate(EXS%MVL,2*NMIDV,Label='EXS%MVL')
      CALL mma_allocate(EXS%MVR,2*NMIDV,Label='EXS%MVR')
      CALL mma_allocate(ISGM,NSGMNT,Label='ISGM')
      CALL mma_allocate(VSGM,NSGMNT,Label='VSGM')
      CALL MKSEG_CP2(SGS%DRT,SGS%DOWN,SGS%LTV,IVR,EXS%MVL,EXS%MVR,ISGM,VSGM,nVert,nLev,NMIDV)

! NIPWLK: NR OF INTEGERS USED TO PACK EACH UP- OR DOWNWALK.
      MXEO=(NLEV*(NLEV+5))/2
      NNRL=(1+MXEO)*NVERT*NSYM
      CALL mma_allocate(EXS%NOCP,MXEO*NMIDV*NSYM,Label='EXS%NOCP')
      CALL mma_allocate(EXS%IOCP,MXEO*NMIDV*NSYM,Label='EXS%IOCP')
      CALL mma_allocate(NRL,NNRL,Label='NRL')
      CALL NRCOUP_CP2(SGS%DRT,ISGM,CIS%NOW,EXS%NOCP,EXS%IOCP,CIS%NOCSF,NRL,EXS%MVL,EXS%MVR,nVert,nMidV,NICOUP,MxEO)
      CALL mma_deallocate(NRL)

      NILNDW=NWALK
      NVTAB_TMP=20000
      NSCR=7*(NLEV+1)
      CALL mma_allocate(EXS%ICOUP,3*NICOUP,Label='EXS%ICOUP')
      CALL mma_allocate(VTAB_TMP,NVTAB_TMP,Label='VTAB_TMP')
      CALL mma_allocate(ILNDW,NILNDW,Label='ILNDW')
      CALL mma_allocate(SCR,NSCR,Label='SCR')
      CALL mma_allocate(VAL,NLEV+1,Label='VAL')
      NICASE=SIZE(CIS%ICASE)
      CALL MKCOUP(nSym,nLev,SGS%ISm,nVert,SGS%MidLev,nMidV,SGS%MVSta,SGS%MVEnd,     &
     &            MxEO,nICoup,nWalk,nICase,nVTAB_TMP,               &
     &            IVR,SGS%MAW,ISGM,VSGM,CIS%NOW,CIS%IOW,EXS%NOCP,EXS%IOCP,ILNDW,      &
     &            CIS%ICase, EXS%ICOUP,VTAB_TMP,NVTAB,SCR,VAL)

      CALL mma_allocate(EXS%VTAB,NVTAB,Label='EXS%VTAB')
      EXS%VTAB(1:NVTAB)=VTAB_TMP(1:NVTAB)

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
      CALL mma_deallocate(SGS%DAW)
      CALL mma_deallocate(SGS%UP)
      CALL mma_deallocate(SGS%RAW)
      Call mma_deallocate(SGS%LTV)

      End Associate

      END SUBROUTINE GINIT_CP2

