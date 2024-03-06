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
      SUBROUTINE SGINIT_CP2(nSym,iSpin,nActEl,nHole1,nEle3,nRas1T,nRas2T,nRas3T,SGS,CIS,EXS,STSYM)
      use stdalloc, only: mma_allocate, mma_deallocate
      use Struct, only: SGStruct, CIStruct, EXStruct
      use gugx, only: IFRAS
      IMPLICIT None
      Integer nSym,iSpin,nActEl,nHole1,nEle3,nRas1T,nRas2T,nRas3T,STSYM
      Type(SGStruct) SGS
      Type(CIStruct) CIS
      Type(EXStruct) EXS
      Integer, Allocatable:: IVR(:), ISGM(:)
      Real*8, Allocatable:: VSGM(:), VTAB_TMP(:)
      Integer NICOUP, NVTAB_TMP, NICASE, nVTab

      Interface
      SUBROUTINE MKGUGA(STSYM,Skip_MKSGNUM)
      IMPLICIT None
      Integer STSYM
      Logical, Optional:: Skip_MKSGNUM
      End SUBROUTINE MKGUGA
      End Interface

      SGS%nSym=nSym
      SGS%iSpin=iSpin
      SGS%nActEl=nActEl

      Associate ( nLev => SGS%nLev, nWalk => CIS%nWalk,                 &
     &            nVert=> SGS%nVert, nMidV=>CIS%nMidV, MXEO => EXS%MxEO, &
     &            LM1RAS=>SGS%LM1RAS, LM3RAS=>SGS%LM3RAS,               &
     &            LV1RAS=>SGS%LV1RAS, LV3RAS=>SGS%LV3RAS,               &
     &            MVSta =>SGS%MVSta,  MVEnd=>SGS%MVEnd)

      LV1RAS=NRAS1T
      LV3RAS=nRas1T+NRAS2T
      LM1RAS=2*nRas1T-NHOLE1
      LM3RAS=NACTEL-NELE3
      IF ((NRAS1T+NRAS3T)/=0) Then
         IFRAS=1
      Else
         IFRAS=0
      End If

      Call mkIsm(SGS)

      Call mknVert0(SGS)

      CALL MKGUGA(STSYM,Skip_MKSGNUM=.TRUE.)

! DECIDE MIDLEV AND CALCULATE MODIFIED ARC WEIGHT TABLE.

      CALL MKMAW(SGS)

! THE DAW, UP AND RAW TABLES WILL NOT BE NEEDED ANY MORE:

! CALCULATE SEGMENT VALUES. ALSO, MVL AND MVR TABLES.
      CALL mma_allocate(IVR,2*NVERT,Label='IVR')
      CALL mma_allocate(EXS%MVL,NMIDV,2,Label='EXS%MVL')
      CALL mma_allocate(EXS%MVR,NMIDV,2,Label='EXS%MVR')
      CALL mma_allocate(ISGM,26*nVert,Label='ISGM')
      CALL mma_allocate(VSGM,26*nVert,Label='VSGM')
      CALL MKSEG(SGS,nVert,nMidV,IVR,EXS%MVL,EXS%MVR,ISGM,VSGM)

! NIPWLK: NR OF INTEGERS USED TO PACK EACH UP- OR DOWNWALK.
      MXEO=(NLEV*(NLEV+5))/2
      CALL mma_allocate(EXS%NOCP,MXEO,NSYM,NMIDV,Label='EXS%NOCP')
      CALL mma_allocate(EXS%IOCP,MXEO,NSYM,NMIDV,Label='EXS%IOCP')
      CALL NRCOUP(SGS,CIS,EXS,nVert,nMidV,MxEO,SGS%ISM,SGS%DRT,ISGM,CIS%NOW,CIS%IOW,EXS%NOCP,EXS%IOCP, &
                      CIS%NOCSF,CIS%IOCSF,CIS%NCSF,EXS%MVL,EXS%MVR,NICOUP,NSYM)

      NVTAB_TMP=20000
      CALL mma_allocate(EXS%ICOUP,3,NICOUP,Label='EXS%ICOUP')
      CALL mma_allocate(VTAB_TMP,NVTAB_TMP,Label='VTAB_TMP')
      NICASE=SIZE(CIS%ICASE)
      CALL MKCOUP(nSym,nLev,SGS%ISm,nVert,SGS%MidLev,nMidV,MVSta,MVEnd,     &
     &            MxEO,nICoup,nWalk,nICase,nVTAB_TMP,               &
     &            IVR,SGS%MAW,ISGM,VSGM,CIS%NOW,CIS%IOW,EXS%NOCP,EXS%IOCP,     &
     &            CIS%ICase, EXS%ICOUP,VTAB_TMP,NVTAB)

      CALL mma_allocate(EXS%VTAB,NVTAB,Label='EXS%VTAB')
      EXS%VTAB(1:NVTAB)=VTAB_TMP(1:NVTAB)

      Call mma_deallocate(VTAB_TMP)
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

      END SUBROUTINE SGINIT_CP2

