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
      Subroutine SGInit(nSym,nActEl,iSpin,nRasPrt,nRas,nRasEl,SGS)
      use stdalloc, only: mma_allocate, mma_deallocate
      use Struct, only: SGStruct
      IMPLICIT None
#include "rassi.fh"
      Integer nSym, nActEl, iSpin, nRasPrt
      Type (SGStruct), Target :: SGS
      ! nRas(iSym,iRasPrt): number of orbitals of symmetry iSym in the iRasPrt active space.
      Integer nRas(8,nRasPrt),nRasEl(nRasPrt)

      Integer, Allocatable:: Lim(:)
      Integer iRO,iSy,Lev

      Interface
      Subroutine MKDRT0(SGS)
      use struct, only: SGStruct
      IMPLICIT None
      Type(SGStruct), Target:: SGS
      End Subroutine MKDRT0
      End Interface

      SGS%nSym=nSym
      SGS%iSpin=iSpin
      SGS%nActEl=nActEl

      Associate ( nLev => SGS%nLev, nVert => SGS%nVert, MidLev => SGS%MidLev, &
                  MVSta => SGS%MVSta, MvEnd => SGS%MVEnd,    &
                  IA0=>SGS%IA0, IB0=>SGS%IB0, IC0=>SGS%IC0, &
                  nVert0=>SGS%nVert0)

      Call mkISM(SGS)

      Call mknVert0(SGS)

! Compute unrestricted DRT tables:
      CALL mma_allocate(SGS%DRT0,NVERT0,5,Label='DRT0')
      CALL mma_allocate(SGS%DOWN0,[1,NVERT0],[0,3],Label='DOWN0')
      nVert=nVert0

      SGS%DRTP => SGS%DRT0
      SGS%DOWNP => SGS%DOWN0

      CALL mkDRT0(SGS)

! Construct a restricted graph.
      Call mma_allocate(Lim,nLev,Label='Lim')
      Lim(:)=0
! Fill in the occupation limit table:
      Lev=0
      Do iRO=1,nRasPrt
        Do iSy=1,nSym
          Lev=Lev+nRas(iSy,iRO)
        End Do
        if(Lev.gt.0) Lim(Lev)=nRasEl(iRO)
      End Do

      Call RmVert(SGS,nLev,Lim)

      Call mma_deallocate(Lim)

      Call mkDRT(SGS)

! Direct Arc Weights table:
      CALL MKDAW(SGS)

! Upchain Index table:
! Reverse Arc Weights table:
      Call MKRAW(SGS)

! Level-To-Vertex table:
      CALL MKLTV(SGS)

      Call MKMID(SGS)

! Modified Arc Weights table:
      CALL MKMAW(SGS)

! The DAW, RAW tables are no longer needed:
      CALL mma_deallocate(SGS%RAW)
      CALL mma_deallocate(SGS%DAW)

      End Associate

      end Subroutine SGInit
