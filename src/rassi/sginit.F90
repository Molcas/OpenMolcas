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
      Subroutine SGInit(nSym,nActEl,iSpin,SGS,CIS)
      use stdalloc, only: mma_deallocate
      use gugx, only: SGStruct, CIStruct
      use MkGUGA_mod, only: MKGUGA
      IMPLICIT None
      Integer nSym, nActEl, iSpin
      Type (SGStruct), Target :: SGS
      Type (CIStruct) :: CIS

      SGS%nSym=nSym
      SGS%iSpin=iSpin
      SGS%nActEl=nActEl

      Call MkGuga(SGS,CIS)

! Modified Arc Weights table:
      CALL MKMAW(SGS)

! The DAW, RAW tables are no longer needed:
      CALL mma_deallocate(SGS%RAW)
      CALL mma_deallocate(SGS%DAW)

      end Subroutine SGInit
