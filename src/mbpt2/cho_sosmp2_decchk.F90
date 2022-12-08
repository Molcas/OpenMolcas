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
! Copyright (C) 2007, Francesco Aquilante                              *
!***********************************************************************

subroutine Cho_SOSmp2_DecChk(irc,iSym,Col,nDim,nCol,Wrk,lWrk,ErrStat)
! Francesco Aquilante, May 2007.
!
! Purpose: check Scaled Opposite-Spin MP2 decomposition of
!          M(ai,bj)=(ai|bj)^2 for sym. block iSym.
!          The columns of the M matrix (actually, the sqrt of
!          their elements) are compared nCol columns at a time
!          with the corresponding integrals obtained from the
!          Cholesky (or RI) vectors.
!          The memory requirement of this routine
!          should be limited to approximately that of the
!          decomposition itself. Note, however, that since all
!          elements are computed, this routine will consume
!          significantly more CPU time.
!          Files are assumed open.
!          On exit,
!          ErrStat(1) = min error
!          ErrStat(2) = max error
!          ErrStat(3) = rms error

use ChoMP2, only: OldVec
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(out) :: irc
integer(kind=iwp), intent(in) :: iSym, nDim, nCol, lWrk
real(kind=wp), intent(inout) :: Col(nDim,nCol)
real(kind=wp), intent(out) :: Wrk(lWrk), ErrStat(3)
integer(kind=iwp) :: iBatCol, ibj1, kai, kbj, lU, Nai, nBatCol, Nbj, NumCol, NumVec
real(kind=wp) :: Fac, xdim
character(len=*), parameter :: SecNam = 'Cho_SOSmp2_DecChk'
real(kind=wp), external :: ddot_
#include "cholesky.fh"
#include "chomp2.fh"
#include "chomp2_dec.fh"

irc = 0

! Check dimensions.
! -----------------

if ((nDim < 1) .or. (nCol < 1)) return
if (nDim /= nT1am(iSym)) then
  irc = -1
  return
end if

! Initialize.
! -----------

ErrStat(1) = huge(ErrStat)
ErrStat(2) = -huge(ErrStat)
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

  ! Compute  M(ai,bj)=(ai|bj)^2  from "new" vectors.
  ! -----------------------------------------------

  lU = lUnit_F(iSym,2)
  NumVec = nMP2Vec(iSym)
  Fac = Zero
  call ChoMP2_DecChk_Int(irc,lU,Col,Nai,Nbj,ibj1,NumVec,Wrk,lWrk,Fac)
  if (irc /= 0) then
    write(u6,*) SecNam,': Cho_SOSmp2_DecChk_Int  rc= ',irc,' [1]'
    irc = 1
    return
  end if

  ! Obtain the (ai|bj) integrals from M(ai,bj).
  ! -------------------------------------------

  do kbj=1,Nbj
    do kai=1,Nai
      Col(kai,kbj) = sqrt(Col(kai,kbj))
    end do
  end do

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
      write(u6,*) SecNam,': Cho_SOSmp2_DecChk_Int returned ',irc,' [2]'
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

xdim = real(Nai*Nai,kind=wp)
ErrStat(3) = sqrt(ErrStat(3)/xdim)

end subroutine Cho_SOSmp2_DecChk
