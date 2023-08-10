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

subroutine Cho_P_SetGL()
!
! Purpose: set global and local index arrays and diagonal.
!          If a sequencial run:
!             Diag   => Diag_Hidden
!             Diag_G => Null()
!          If a parallel run:
!             Diag   => Diag_G_Hidden
!             Diag_G => Diag_Hidden
!
!          Diag_Hidden is allocated in the calling routine, while
!          Diag_G_Hidden is allocated here if needed.

use Cholesky, only: Cho_Real_Par, Diag, Diag_G, Diag_G_Hidden, Diag_Hidden, iiBstR, iiBstR_G, iiBstRSh, iiBstRSh_G, &
                    iiBstRsh_L_Hidden, iL2G, IndRed, IndRed_G, IndRed_G_Hidden, IndRSh, IndRSh_G, IndRsh_G_Hidden, InfRed, &
                    InfRed_G, InfRed_G_Hidden, InfVec, InfVec_G, InfVec_G_Hidden, LuPri, mmBstRT, mmBstRT_G, MySP, n_MySP, nnBstR, &
                    nnBstR_G, nnBstRSh, nnBstRSh_G, nnBstRT, nnBstRT_G, nnBstRsh_L_Hidden, nnShl, nnShl_G, nSym
use stdalloc, only: mma_allocate
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: i, i1, i2, irc, iShlAB, iSP, iSym, N
character(len=*), parameter :: SecNam = 'Cho_P_SetGL'

! If not parallel, return.
! ------------------------

if (.not. Cho_Real_Par) then
  Diag => Diag_Hidden
  return
end if

! Set global data Cholesky.
! -------------------------

Diag_G => Diag_Hidden

nnShl_G = nnShl
mmBstRT_G = mmBstRT

iiBstR_G(:,:) = iiBstR(:,:)
nnBstR_G(:,:) = nnBstR(:,:)
nnBstRT_G(:) = nnBstRT(:)

InfRed_G => InfRed

InfVec_G => InfVec

iiBstRSh_G => iiBstRSh

nnBstRSh_G => nnBstRSh

IndRed_G => IndRed

IndRSh_G => IndRSh

! Reallocate and reset local data.
! --------------------------------

call mma_allocate(InfRed_G_Hidden,size(InfRed),Label='InfRed_G_Hidden')
InfRed => InfRed_G_Hidden

call mma_allocate(InfVec_G_Hidden,size(InfVec,1),size(InfVec,2),size(InfVec,3),Label='InfVec_G_Hidden')
InfVec => InfVec_G_Hidden

nnShl = n_mySP
call mma_allocate(iiBstRsh_L_Hidden,nSym,n_mySP,3,Label='iiBstRSh_L_Hidden')
iiBstRSh => iiBstRSh_L_Hidden
call mma_allocate(nnBstRsh_L_Hidden,nSym,n_mySP,3,Label='nnBstRSh_L_Hidden')
nnBstRSh => nnBstRSh_L_Hidden

do iSP=1,nnShl
  iShlAB = mySP(iSP)
  nnBstRSh(1:nSym,iSP,1) = nnBstRSh_G(1:nSym,iShlAB,1)
end do
call Cho_SetRedInd(1)
mmBstRT = nnBstRT(1)

call mma_allocate(IndRed_G_Hidden,mmBstRT,3,Label='IndRed_G_Hidden')
IndRed => IndRed_G_Hidden
call mma_allocate(IndRSh_G_Hidden,mmBstRT,Label='IndRSh_G_Hidden')
IndRSh => IndRSh_G_Hidden
call mma_allocate(iL2G,mmBstRT,Label='iL2G')

N = 0
do iSym=1,nSym
  do iSP=1,nnShl
    iShlAB = mySP(iSP)
    i1 = iiBstR_G(iSym,1)+iiBstRSh_G(iSym,iShlAB,1)+1
    i2 = i1+nnBstRSh_G(iSym,iShlAB,1)-1
    do i=i1,i2
      N = N+1
      IndRed(N,1) = IndRed_G(i,1)
      IndRSh(N) = IndRSh_G(i)
      iL2G(N) = i
    end do
  end do
end do
call Cho_X_RSCopy(irc,1,2)
if (irc /= 0) then
  write(Lupri,*) SecNam,': [1] Cho_X_RSCopy returned ',irc
  call Cho_Quit('Error in '//SecNam,104)
end if
call Cho_X_RSCopy(irc,2,3)
if (irc /= 0) then
  write(Lupri,*) SecNam,': [2] Cho_X_RSCopy returned ',irc
  call Cho_Quit('Error in '//SecNam,104)
end if

! Allocate and set local diagonal.
! --------------------------------

call mma_allocate(Diag_G_Hidden,mmBstRT,Label='Diag_G_Hidden')
Diag => Diag_G_Hidden
do i=1,mmBstRT
  Diag(i) = Diag_G(iL2G(i))
end do

end subroutine Cho_P_SetGL
