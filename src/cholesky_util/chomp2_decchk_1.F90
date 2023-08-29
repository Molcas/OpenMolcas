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
! Copyright (C) 2005, Thomas Bondo Pedersen                            *
!***********************************************************************

subroutine ChoMP2_DecChk_1(irc,iSym,Col,nDim,nCol,Wrk,lWrk,ErrStat)
!
! Thomas Bondo Pedersen, Jan. 2005.
!
! Purpose: check MP2 decomposition of the (ai|bj) integrals
!          (sym. block iSym). The columns of the (ai|bj) matrix are
!          compared nCol columns at a time. This
!          implies that the memory requirement of this routine
!          should be limited to approximately that of the
!          decomposition itself. Note, however, that since all
!          integrals are computed, this routine will consume
!          significantly more CPU time.
!          Files are assumed open.
!          On exit,
!          ErrStat(1) = min error
!          ErrStat(2) = max error
!          ErrStat(3) = rms error

use Cholesky, only: NumCho
use ChoMP2, only: Incore, lUnit_F, nMP2Vec, nT1am, OldVec
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(out) :: irc
integer(kind=iwp), intent(in) :: iSym, nDim, nCol, lWrk
real(kind=wp), intent(inout) :: Col(nDim,nCol)
real(kind=wp), intent(out) :: Wrk(lWrk), ErrStat(3)
integer(kind=iwp) :: iBatCol, ibj1, kai, kbj, lU, Nai, nBatCol, Nbj, NumCol, NumVec
real(kind=wp) :: Fac, xdim
character(len=*), parameter :: SecNam = 'ChoMP2_DecChk_1'
real(kind=wp), external :: ddot_

irc = 0

! Check dimensions.
! -----------------

if ((nDim < 1) .or. (nCol < 1)) return
if (nDim /= nT1am(iSym)) then
  irc = -1
  return ! exit
end if

! Initialize.
! -----------

ErrStat(1) = 9.9e15_wp
ErrStat(2) = -9.9e15_wp
ErrStat(3) = Zero

! Set up batching over columns of the (ai|bj) matrix.
! ---------------------------------------------------

NumCol = min(nCol,nT1am(iSym))
nBatCol = (nT1am(iSym)-1)/NumCol+1

! Start column batch loop.
! ------------------------

Nai = nDim
do iBatCol=1,nBatCol

  ! Set batch info.
  ! ---------------

  if (iBatCol == nBatCol) then
    Nbj = nT1am(iSym)-NumCol*(nBatCol-1)
  else
    Nbj = NumCol
  end if
  ibj1 = NumCol*(iBatCol-1)+1

  ! Compute integrals from "new" vectors.
  ! -------------------------------------

  lU = lUnit_F(iSym,2)
  NumVec = nMP2Vec(iSym)
  Fac = Zero
  call ChoMP2_DecChk_Int(irc,lU,Col,Nai,Nbj,ibj1,NumVec,Wrk,lWrk,Fac)
  if (irc /= 0) then
    write(u6,*) SecNam,': ChoMP2_DecChk_Int returned ',irc,' [1]'
    irc = 1
    return
  end if

  ! Compute "old" and subtract "new".
  ! ---------------------------------

  if (InCore(iSym)) then
    call DGEMM_('N','T',Nai,Nbj,NumCho(iSym),One,OldVec,Nai,OldVec(ibj1),Nai,-One,Col,Nai)
  else
    lU = lUnit_F(iSym,1)
    NumVec = NumCho(iSym)
    Fac = -One
    call ChoMP2_DecChk_Int(irc,lU,Col,Nai,Nbj,ibj1,NumVec,Wrk,lWrk,Fac)
    if (irc /= 0) then
      write(u6,*) SecNam,': ChoMP2_DecChk_Int returned ',irc,' [2]'
      irc = 2
      return
    end if
  end if

  ! Compute error stats.
  ! --------------------

  do kbj=1,Nbj
    do kai=1,Nai
      ErrStat(1) = min(ErrStat(1),Col(kai,kbj))
      ErrStat(2) = max(ErrStat(2),Col(kai,kbj))
    end do
  end do
  ErrStat(3) = ErrStat(3)+dDot_(Nai*Nbj,Col,1,Col,1)

end do

! Compute rms error.
! ------------------

xdim = real(Nai,kind=wp)**2
ErrStat(3) = sqrt(ErrStat(3)/xdim)

end subroutine ChoMP2_DecChk_1
