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

      Real*8,  Allocatable:: VTabTmp(:)
      Integer nVMax, nVTab, nVTab_final

      Associate ( nLev => SGS%nLev, nVert => SGS%nVert,                 &
     &            MidLev =>SGS%MidLev, MVSta => SGS%MVSta,              &
     &            MVEnd  =>SGS%MVEnd, nMidV=>CIS%nMidV,                 &
     &            nIpWlk => CIS%nIpWlk, MxEO=>EXS%MxEO,                 &
     &            nWalk=>CIS%nWalk)

      nMidV=MVEnd-MVSta+1
      nIpWlk=1+(MidLev-1)/15
      nIpWlk=MAX(nIpWlk,1+(nLev-MidLev-1)/15)

! Calculate segment values, and MVL and MVR tables:

      Call MkSeg(SGS,CIS,EXS)


! Various offset tables:

      Call NrCoup(SGS,CIS,EXS)

! Computed in NrCoup:

      nVMax=5000
      Call mma_allocate(VTabTmp,nVMax,Label='VTabTmp')
      nVTab=nVMax
      Call MkCoup(MidLev,MVSta,MVEnd,nWalk,nVTab,                       &
     &            VTabTmp,NVTAB_Final,SGS,CIS,EXS)

      nVTAB=nVTAB_Final
      Call mma_allocate(EXS%Vtab,nVTab,Label='EXS%VTab')
      EXS%VTab(1:nVTab)=VTabTmp(1:nVTab)
      Call mma_deallocate(VTabTmp)

      Call mma_deallocate(CIS%ISgm)
      Call mma_deallocate(CIS%VSgm)
      Call mma_deallocate(CIS%IVR)

      End Associate

      end Subroutine CXInit
