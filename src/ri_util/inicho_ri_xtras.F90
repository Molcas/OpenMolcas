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

subroutine IniCho_RI_Xtras(iTOffs,nIrrep,iShij,nShij)

use Cholesky, only: iiBstR, iiBstRSh, IndRed, IndRed_Hidden, IndRSh, IndRSh_Hidden, iRS2F, MaxRed, mmBstRT, nnBstR, nnBstRSh, &
                    nnBstRT, nnShl, nDimRS, nSym
use stdalloc, only: mma_allocate
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: nIrrep, iTOffs(3,nIrrep), nShij, iShij(2,nShij)
integer(kind=iwp) :: i, iCount, iiBst(8), iSym, nnBst(8), nnBstT
logical(kind=iwp) :: DoDummy

! Define max. dimensions and offsets of the symmetry blocks of the
! integrals matrix.
! ----------------------------------------------------------------

iCount = 0
do iSym=1,nSym
  iiBst(iSym) = iCount
  nnBst(iSym) = iTOffs(3,iSym)
  iCount = iCount+nnBst(iSym)
end do

! Set dimensions of reduced sets equal to full dimension.
! -------------------------------------------------------

do i=1,3
  nnBstT = 0
  do iSym=1,nSym
    iiBstR(iSym,i) = iiBst(iSym)
    nnBstR(iSym,i) = nnBst(iSym)
    nnBstT = nnBstT+nnBstR(iSym,i)
  end do
  nnBstRT(i) = nnBstT
end do
mmBstRT = nnBstRT(1)

! Allocate index arrays for reduced sets: IndRed and IndRsh.
! ----------------------------------------------------------

call mma_allocate(IndRed_Hidden,nnBstRT(1),3,Label='IndRed_Hidden')
IndRed => IndRed_Hidden
call mma_allocate(IndRSh_Hidden,nnBstRT(1),Label='IndRSh_Hidden')
IndRSh => IndRSh_Hidden

! Allocate iScr array used by reading routines.
! ---------------------------------------------

DoDummy = .false.
call Cho_Allo_iScr(DoDummy)

! Initialize reduced set dimensions used for reading vectors.
! (Note: here they are all the same - there is one reduced sets!)
! ---------------------------------------------------------------

do i=1,MaxRed
  nDimRS(1:nSym,i) = nnBstR(1:nSym,1)
end do

! Allocate and set mapping array from 1st reduced set to full storage.
! --------------------------------------------------------------------

call mma_allocate(iRS2F,2,nnBstRT(1),Label='iRS2F')

! Set index arrays corresponding to full storage:
! iiBstRSh, nnBstRSh, IndRed, IndRSh, and iRS2F.
! -----------------------------------------------

call SetChoIndx_RI(iiBstRSh,nnBstRSh,IndRed,IndRsh,iRS2F,nSym,nnShl,nnBstRT(1),iShij,nShij)

return

end subroutine IniCho_RI_Xtras
