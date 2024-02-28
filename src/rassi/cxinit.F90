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
      Subroutine CXInit(SGS,CIS,EXS,nSym)
      use stdalloc, only: mma_allocate, mma_deallocate
      use Struct, only:  SGStruct, CIStruct, EXStruct
      IMPLICIT None
      Integer nSym
      Type (SGStruct) SGS
      Type (CIStruct) CIS
      Type (EXStruct) EXS

      Integer, Allocatable:: IVR(:), ISgm(:), NRL(:), iLndw(:), Scr(:)
      Real*8,  Allocatable::         VSgm(:)
      Real*8,  Allocatable:: VTabTmp(:), Val(:)
      Integer nLev, nVert, MidLev, MVSta, MVEnd, nMidV, nIpWlk, nSgmnt, &
     &        MxEO, nNRL, nWalk,    &
     &        nICoup, nVMax, niLndw, nICase, nScr, nVTab,      &
     &        nVTab_final

      nLev   =SGS%nLev
      nVert  =SGS%nVert
      MidLev =SGS%MidLev
      MVSta  =SGS%MVSta
      MVEnd  =SGS%MVEnd

! Calculate segment values, and MVL and MVR tables:
      nMidV=1+MVEnd-MVSta
! nIpWlk: NR OF INTEGERS USED TO PACK EACH UP- OR DOWNWALK.
      nIpWlk=1+(MidLev-1)/15
      nIpWlk=MAX(nIpWlk,1+(nLev-MidLev-1)/15)
      Call mma_allocate(IVR,2*nVert,Label='IVR')
      Call mma_allocate(EXS%MVR,nMidV,2,Label='EXS%MVR')
      Call mma_allocate(EXS%MVL,nMidV,2,Label='EXS%MVL')
      nSgmnt=26*nVert
      Call mma_allocate(ISgm,nSgmnt,Label='ISgm')
      Call mma_allocate(VSgm,nSgmnt,Label='VSgm')
      Call MkSeg(SGS,nLev,nVert,nMidv,SGS%DRT,SGS%Down,SGS%LTV,IVR,EXS%MVL,EXS%MVR,ISgm,VSgm)
      CIS%nMidV   =nMidV
      CIS%nIpWlk  = nIpWlk

! Various offset tables:
      Call mma_allocate(CIS%NOW,2,nSym,nMidV,Label='CIS%NOW')
      Call mma_allocate(CIS%IOW,2,nSym,nMidV,Label='CIS%IOW')
      MxEO=(nLev*(nLev+5))/2
      nNRL=(1+MxEO)*nVert*nSym
      Call mma_allocate(EXS%NOCP,MxEO,nSym,nMidV,Label='EXS%NOCP')
      Call mma_allocate(EXS%IOCP,MxEO,nSym,nMidV,Label='EXS%IOCP')
      Call mma_allocate(CIS%NCSF,nSym,Label='CIS%NCSF')

      Call mma_allocate(CIS%NOCSF,nSym,nMidV,nSym,Label='CIS%NOCSF')
      Call mma_allocate(CIS%IOCSF,nSym,nMidV,nSym,Label='CIS%IOCSF')
      EXS%MxEO =MxEO
      Call mma_allocate(NRL,nNRL,Label='NRL')
      Call NrCoup(SGS,CIS,EXS,nVert,nMidV,MxEO,SGS%ISm,SGS%DRT,         &
     &            ISgm,CIS%NOW,CIS%IOW,EXS%NOCP,                        &
     &            EXS%IOCP,CIS%NOCSF,CIS%IOCSF,                         &
     &            CIS%NCSF,NRL,EXS%MVL,EXS%MVR,nICoup,nSym)
      Call mma_deallocate(NRL)
! Computed in NrCoup:
      nWalk=CIS%nWalk

      nICase=nWalk*nIpWlk
      Call mma_allocate(CIS%ICase,nICase,Label='CIS%ICase')
      Call mma_allocate(EXS%ICoup,3,nICoup,Label='EXS%ICoup')
      nVMax=5000
      Call mma_allocate(VTabTmp,nVMax,Label='VTabTmp')
      nILNDW=nWalk
      Call mma_allocate(iLndw,niLndw,Label='iLndw')
      nScr=7*(nLev+1)
      Call mma_allocate(Scr,nScr,Label='Scr')
      Call mma_allocate(Val,nLev+1,Label='Val')
      nVTab=nVMax
      Call MkCoup(nSym,nLev,SGS%Ism,nVert,MidLev,nMidV,MVSta,MVEnd,     &
     &            MxEO,nICoup,nWalk,nICase,nVTab,                       &
     &            IVR,SGS%MAW,ISGM,VSGM,CIS%NOW,CIS%IOW,EXS%NOCP,       &
     &            EXS%IOCP,ILNDW,CIS%ICase,EXS%ICOUP,VTabTmp,           &
     &            NVTAB_Final,SCR,VAL)

      nVTAB=nVTAB_Final
      Call mma_allocate(EXS%Vtab,nVTab,Label='EXS%VTab')
      EXS%VTab(1:nVTab)=VTabTmp(1:nVTab)
      Call mma_deallocate(VTabTmp)

      Call mma_deallocate(iLndw)
      Call mma_deallocate(Scr)
      Call mma_deallocate(Val)
      Call mma_deallocate(ISgm)
      Call mma_deallocate(VSgm)
      Call mma_deallocate(IVR)

      end Subroutine CXInit
