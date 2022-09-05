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
! Copyright (C) Thomas Bondo Pedersen                                  *
!***********************************************************************
      SubRoutine IniCho_RI_Xtras(iTOffs,nIrrep,iShij,nShij)
      use ChoArr, only: iRS2F, nDimRS
      use ChoSwp, only: nnBstRSh, iiBstRSh
      use ChoSwp, only:   IndRSh,   IndRSh_Hidden
      use ChoSwp, only:   IndRed,   IndRed_Hidden
      Implicit None
      Integer nIrrep, nShij
#include "cholesky.fh"
#include "choorb.fh"
#include "stdalloc.fh"

      Logical DoDummy

      Integer iiBst(8), nnBst(8), iTOffs(3,nIrrep), iShij(2,nShij)
      Integer iSym, iCount, nnBstT

      Integer i

!     Define max. dimensions and offsets of the symmetry blocks of the
!     integrals matrix.
!     ----------------------------------------------------------------

      iCount = 0
      Do iSym = 1,nSym
         iiBst(iSym) = iCount
         nnBst(iSym) = iTOffs(3,iSym)
         iCount = iCount + nnBst(iSym)
      End Do

!     Set dimensions of reduced sets equal to full dimension.
!     -------------------------------------------------------

      Do i = 1,3
         nnBstT=0
         Do iSym = 1,nSym
            iiBstR(iSym,i) = iiBst(iSym)
            nnBstR(iSym,i) = nnBst(iSym)
            nnBstT=nnBstT+nnBstR(iSym,i)
         End Do
         nnBstRT(i) = nnBstT
      End Do
      mmBstRT = nnBstRT(1)

!     Allocate index arrays for reduced sets: IndRed and IndRsh.
!     ----------------------------------------------------------

      Call mma_allocate(IndRed_Hidden,nnBstRT(1),3,                     &
     &                  Label='IndRed_Hidden')
      IndRed => IndRed_Hidden
      Call mma_allocate(IndRSh_Hidden,nnBstRT(1),Label='IndRSh_Hidden')
      IndRSh => IndRSh_Hidden

!     Allocate iScr array used by reading routines.
!     ---------------------------------------------

      DoDummy = .False.
      Call Cho_Allo_iScr(DoDummy)

!     Initialize reduced set dimensions used for reading vectors.
!     (Note: here they are all the same - there is one reduced sets!)
!     ---------------------------------------------------------------

      Do i = 1,MaxRed
         Do iSym = 1,nSym
            nDimRS(iSym,i) = nnBstR(iSym,1)
         End Do
      End Do

!     Allocate and set mapping array from 1st reduced set to full
!     storage.
!     -----------------------------------------------------------

      Call mma_allocate(iRS2F,2,nnBstRT(1),Label='iRS2F')

!     Set index arrays corresponding to full storage:
!     iiBstRSh, nnBstRSh, IndRed, IndRSh, and iRS2F.
!     -----------------------------------------------

      Call SetChoIndx_RI(iiBstRSh,nnBstRSh,                             &
     &                   IndRed,IndRsh,iRS2F,                           &
     &                   nSym,nnShl,nnBstRT(1),3,2,iShij,nShij)

      Return
      End
