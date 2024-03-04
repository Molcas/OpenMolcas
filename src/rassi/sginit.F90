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
      Subroutine SGInit(nSym,nActEl,iSpin,SGS)
      use stdalloc, only: mma_allocate, mma_deallocate
      use Struct, only: SGStruct
      IMPLICIT None
#include "rassi.fh"
      Integer nSym, nActEl, iSpin
      Type (SGStruct), Target :: SGS

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
                  MVSta => SGS%MVSta, MvEnd => SGS%MVEnd, nVert0=>SGS%nVert0)

      Call mkISM(SGS)

      Call mknVert0(SGS)

! Compute unrestricted DRT tables:
      CALL mma_allocate(SGS%DRT0,NVERT0,5,Label='DRT0')
      CALL mma_allocate(SGS%DOWN0,[1,NVERT0],[0,3],Label='DOWN0')
      nVert=nVert0

      SGS%DRTP => SGS%DRT0
      SGS%DOWNP => SGS%DOWN0

      CALL mkDRT0(SGS)

      Call mkRAS(SGS)

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
