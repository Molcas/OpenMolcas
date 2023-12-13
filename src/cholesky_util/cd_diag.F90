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
! Copyright (C) 2004, Thomas Bondo Pedersen                            *
!***********************************************************************

subroutine CD_Diag(CD_Vec,Restart,Converged,Thr,ThrNeg,ThrFail,DiaInp,Diag,Buf,nDim,lBuf,ErrStat,NumCho,irc)
!
! Thomas Bondo Pedersen, October 2004.
!
! Purpose: set up diagonal for general Cholesky decomposition.
!
! NB!! ThrNeg and ThrFail are supposed to be negative with
!      ThrFail < ThrNeg. Diagonals less than ThrNeg are zeroed,
!      while diagonals less than ThrFail is taken as a sign of
!      a non-positive definite matrix (i.e., decomposition will
!      fail).
!
! Error codes, irc:
!    0 : all OK
!  201 : inconsistent input: Restart but NumCho < 0.
!        (NumCho >= 0 is allowed with Restart.)
!  202 : insufficient buffer size, lBuf.
!  203 : too negative diagonal element found (i.e., matrix
!        is non-positive definite).
!
! CD_Vec : external routine for vectors

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
external :: CD_Vec
logical(kind=iwp), intent(in) :: Restart
logical(kind=iwp), intent(out) :: Converged
integer(kind=iwp), intent(in) :: nDim, lBuf, NumCho
real(kind=wp), intent(in) :: Thr, ThrNeg, ThrFail, DiaInp(nDim)
real(kind=wp), intent(out) :: Diag(nDim), ErrStat(3)
real(kind=wp), intent(inout) :: Buf(lBuf)
integer(kind=iwp), intent(out) :: irc
integer(kind=iwp) :: i, iBatch, ij, iOpt, iVec1, jVec, kOff, nBatch, NumV, nVec
real(kind=wp) :: xDim

! Set variables.
! --------------

irc = 0
if (nDim < 1) then
  Converged = .true. ! in a sense, at least
  return ! exit (nothing to do)
else
  Converged = .false.
end if

! Set up diagonal.
! ----------------

Diag(:) = DiaInp(:)

if (Restart) then ! subtract previous vectors
  if (NumCho > 0) then

    nVec = min(NumCho,lBuf/nDim)
    if (nVec < 1) then
      irc = 202
      return ! exit (insufficient buffer size)
    else
      nBatch = (NumCho-1)/nVec+1
    end if

    do iBatch=1,nBatch

      if (iBatch == nBatch) then
        NumV = NumCho-nVec*(nBatch-1)
      else
        NumV = nVec
      end if

      iVec1 = nVec*(iBatch-1)+1
      iOpt = 2
      call CD_Vec(iVec1,NumV,Buf,lBuf,nDim,iOpt)

      do jVec=1,NumV
        kOff = nDim*(jVec-1)
        do i=1,nDim
          ij = kOff+i
          Diag(i) = Diag(i)-Buf(ij)*Buf(ij)
        end do
      end do

    end do

  else if (NumCho < 0) then

    irc = 201
    return ! exit (inconsistent input)

  end if
end if

! Zero negative diagonals (fail if too negative),
! get error statistics (min, max, and rms error),
! check convergence based on max diagonal.
! -----------------------------------------------

if (Diag(1) < ThrNeg) then
  if (Diag(1) < ThrFail) then
    irc = 203
    return ! exit (too negative diagonal)
  else
    Diag(1) = Zero
  end if
end if
ErrStat(1) = Diag(1)
ErrStat(2) = Diag(1)
ErrStat(3) = Diag(1)*Diag(1)
do i=2,nDim
  if (Diag(1) < ThrNeg) then
    if (Diag(1) < ThrFail) then
      irc = 203
      return ! exit (too negative diagonal)
    else
      Diag(1) = Zero
    end if
  end if
  ErrStat(1) = min(ErrStat(1),Diag(i))
  ErrStat(2) = max(ErrStat(2),Diag(i))
  ErrStat(3) = ErrStat(3)+Diag(i)*Diag(i)
end do
xDim = real(nDim,kind=wp)
ErrStat(3) = sqrt(ErrStat(3))/xDim

Converged = ErrStat(2) <= Thr

end subroutine CD_Diag
