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

subroutine Cho_X_GenVec(irc,Diag)

use Cholesky, only: iQuAB, pTemp, iQuAB_here, LuPri, MaxQual, Mode_Screen, nnZTot, nSym, NumCho
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(out) :: irc
real(kind=wp), intent(_OUT_) :: Diag(*)
integer(kind=iwp) :: iSym, MaxQual_SAVE
character(len=*), parameter :: SecNam = 'Cho_X_GenVec'

! Set return code.
! ----------------

irc = 0

! Re-allocate the iQuAB index array, save old allocation.
! This is used to trick the integral extraction from Seward so as to
! reduce the number of re-calculations of shell pairs.
! ------------------------------------------------------------------

pTemp => iQuAB
MaxQual_SAVE = MaxQual

MaxQual = NumCho(1)
do iSym=2,nSym
  MaxQual = max(MaxQual,NumCho(iSym))
end do

call mma_allocate(iQuAB_here,MaxQual,nSym,Label='iQuAB_here')
iQuAB => iQuAB_here

! Read initial diagonal.
! ----------------------

call Cho_IODiag(Diag,2)

! Reinitialize the number of zeroed negative diagonals.
! Turn on damped screening for second step.
! -----------------------------------------------------

nNZTot = 0
MODE_SCREEN = 1

! Generate vectors.
! -----------------

call Cho_GnVc_Drv(irc,Diag)
if (irc /= 0) write(Lupri,*) SecNam,': Cho_GnVc_Drv returned ',irc

! De-allocations.
! Restore original iQuAB array.
! -----------------------------

call mma_deallocate(iQuAB_here)
iQuAB => pTemp
MaxQual = MaxQual_SAVE

end subroutine Cho_X_GenVec
